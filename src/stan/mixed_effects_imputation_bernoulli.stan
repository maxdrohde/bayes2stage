data {
  // ... standard integers ...
  int<lower=0> N;
  int<lower=1> G;
  int<lower=0, upper=G> G_obs;
  int<lower=0, upper=G> G_mis;

  // ... standard data ...
  int<lower=0> P;
  vector[N] t;
  matrix[N, P] X;
  int<lower=0> S;
  matrix[G, S] Z;
  vector[N] y;

  // *** OPTIMIZATION REQUIREMENTS ***
  // 1. Array to broadcast G random effects to N observations
  array[N] int<lower=1, upper=G> id; 

  // 2. Arrays to sum N observations back to G subjects (Data must be sorted!)
  array[G] int<lower=1> pos; 
  array[G] int<lower=0> len; 

  // ... missing data indices ...
  array[G_obs] int<lower=1, upper=G> index_obs;
  array[G_mis] int<lower=1, upper=G> index_mis;
  array[G_obs] int<lower=0, upper=1> x_obs;
}

parameters {
  // ... (parameters remain exactly the same) ...
  matrix[2, G] z_re; 
  vector<lower=0>[2] sigma_re;
  cholesky_factor_corr[2] L_re; 

  real alpha_main;
  real beta_t;
  real beta_x;
  real beta_x_t_interaction;
  vector[P] beta;
  real<lower=0> sigma_main;

  real alpha_imputation;
  vector[S] gamma;
}

model {
  // --- Priors ---
  to_vector(z_re) ~ std_normal();
  sigma_re ~ exponential(0.1);
  L_re ~ lkj_corr_cholesky(2);

  sigma_main ~ exponential(0.1);
  alpha_main ~ normal(0, 100);
  beta ~ normal(0, 100);
  beta_t ~ normal(0, 100);
  beta_x ~ normal(0, 100);
  beta_x_t_interaction ~ normal(0, 100);

  alpha_imputation ~ normal(0, 100);
  gamma ~ normal(0, 100);

  // --- 1. Linear Algebra (Vectorized over N) ---
  
  // Construct Random Effects Matrix (2 x G)
  matrix[2, G] re = diag_pre_multiply(sigma_re, L_re) * z_re;

  // BROADCASTING: Expand Random Effects from G to N using 'id'
  // This is extremely fast in Stan
  vector[N] re_intercept_N = re[1, id]'; 
  vector[N] re_slope_N     = re[2, id]';

  // Fixed Effects (N)
  vector[N] xb = X * beta;

  // Base Linear Predictor (N)
  // This calculates the mean for every observation in one giant operation
  vector[N] base = alpha_main 
                   + re_intercept_N 
                   + (beta_t + re_slope_N) .* t 
                   + xb;

  // Interaction term (N)
  vector[N] delta = beta_x + beta_x_t_interaction * t;

  // --- 2. Manual Log-Likelihood Calculation (Vectorized) ---
  // Instead of calling normal_lpdf G times, we calculate the densities manually 
  // for all N. This avoids function call overhead.
  // Log Normal: -0.5 * log(2pi) - log(sigma) - 0.5 * ((y - mu)/sigma)^2

  real inv_sigma2 = inv(square(sigma_main));
  real log_sigma  = log(sigma_main);
  real const_term = -0.5 * log(2 * pi());

  // Calculate squared residuals for the two cases (x=0 and x=1)
  // These are vectors of length N
  vector[N] res_sq_x0 = square(y - base);
  vector[N] res_sq_x1 = square(y - (base + delta));

  // --- 3. Aggregation (Summing N back to G) ---
  
  vector[G] eta_imputation = alpha_imputation + Z * gamma;
  vector[G] log_y_x0;
  vector[G] log_y_x1;

  // We loop G times only to sum the pre-calculated residuals.
  // This is very cheap because no new graph nodes are created here, just summation.
  for (g in 1:G) {
    int s = pos[g];
    int l = len[g];
    int e = s + l - 1;

    // Sum of squared residuals for this subject
    real sum_sq_0 = sum(res_sq_x0[s:e]);
    real sum_sq_1 = sum(res_sq_x1[s:e]);

    // Construct the final log-likelihood for this subject
    // LL = N_g * (const - log_sigma) - 0.5 * (sum_sq / sigma^2)
    real ll_const = l * (const_term - log_sigma);

    log_y_x0[g] = ll_const - 0.5 * sum_sq_0 * inv_sigma2;
    log_y_x1[g] = ll_const - 0.5 * sum_sq_1 * inv_sigma2;
  }

  // --- 4. Joint Likelihood (Same as before) ---

  // Observed X
  for (k in 1:G_obs) {
    int g = index_obs[k];
    int xg = x_obs[k];
    target += bernoulli_logit_lpmf(xg | eta_imputation[g]);
    target += (xg == 1) ? log_y_x1[g] : log_y_x0[g];
  }

  // Missing X (Marginalization)
  for (k in 1:G_mis) {
    int g = index_mis[k];
    target += log_sum_exp(
      bernoulli_logit_lpmf(1 | eta_imputation[g]) + log_y_x1[g],
      bernoulli_logit_lpmf(0 | eta_imputation[g]) + log_y_x0[g]
    );
  }
}

generated quantities {
  corr_matrix[2] corr_rand_effects = multiply_lower_tri_self_transpose(L_re);
  matrix[2, 2] cov_rand_effects = quad_form_diag(corr_rand_effects, sigma_re);
}