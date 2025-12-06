data {
  int<lower=0> N;                   // number of observations
  int<lower=1> G;                   // number of subjects
  int<lower=0, upper=G> G_obs;      // number of subjects with observed x
  int<lower=0, upper=G> G_mis;      // number of subjects with missing x

  int<lower=0> P;                   // number of main model predictors
  vector[N] t;                      // time
  matrix[N, P] X;                   // main model covariates

  int<lower=0> S;                   // number of imputation model predictors
  matrix[G, S] Z;                   // imputation model covariates

  vector[N] y;                      // outcome

  array[G_obs] int<lower=1, upper=G> index_obs;
  array[G_mis] int<lower=1, upper=G> index_mis;

  // Count covariate for subjects with observed x
  array[G_obs] int<lower=0> x_obs;
  
  // NEW: A fixed upper bound for marginalization (e.g., 200, 500)
  // Should be large enough that P(x > max_x) is negligible.
  int<lower=max(x_obs)> max_x;

  // Subject index per row
  array[N] int<lower=1, upper=G> id;
}

parameters {
  // Non-centered random effects
  matrix[2, G] z_re;
  vector<lower=0>[2] sigma_re;
  cholesky_factor_corr[2] L_re;

  // Main model
  real alpha_main;
  real beta_t;
  real beta_x;                       
  real beta_x_t_interaction;
  vector[P] beta;
  real<lower=0> sigma_main;

  // Imputation model for negative binomial x
  real alpha_imputation;
  vector[S] gamma;
  real<lower=0> phi;                 
}

model {
  // === Priors ===
  // Random effects
  to_vector(z_re) ~ std_normal();
  sigma_re ~ exponential(0.1);
  L_re ~ lkj_corr_cholesky(2);

  // Main model
  sigma_main ~ exponential(0.1);
  alpha_main ~ normal(0, 2);
  beta ~ normal(0, 2);
  beta_t ~ normal(0, 2);
  beta_x ~ normal(0, 2);
  beta_x_t_interaction ~ normal(0, 2);

  // Imputation model
  alpha_imputation ~ normal(0, 2.5);
  gamma ~ normal(0, 2.5);
  phi ~ exponential(0.1);

  // === Transformed quantities (Local) ===
  matrix[2, G] re = diag_pre_multiply(sigma_re, L_re) * z_re;
  vector[G] mu_imp = exp(alpha_imputation + Z * gamma);
  vector[N] xb = X * beta;

  // === Precompute sufficient statistics ===
  vector[G] N_g = rep_vector(0, G);   
  vector[G] A_g = rep_vector(0, G);   
  vector[G] B_g = rep_vector(0, G);   
  vector[G] C_g = rep_vector(0, G);   

  for (n in 1:N) {
    int g = id[n];
    real base = alpha_main + re[1, g] + (beta_t + re[2, g]) * t[n] + xb[n];
    real resid = y[n] - base;
    real delta = beta_x + beta_x_t_interaction * t[n];

    N_g[g] += 1;
    A_g[g] += square(resid);
    B_g[g] += resid * delta;
    C_g[g] += square(delta);
  }

  real log_2pi_sigma2 = log(2 * pi() * square(sigma_main));
  real inv_2sigma2 = 0.5 / square(sigma_main);

  // === Likelihood: Subjects with Observed X ===
  for (k in 1:G_obs) {
    int g = index_obs[k];
    int xg = x_obs[k];
    real x_val = xg;

    // 1. Imputation model
    target += neg_binomial_2_lpmf(xg | mu_imp[g], phi);

    // 2. Main model
    real SS = A_g[g] - 2 * B_g[g] * x_val + C_g[g] * square(x_val);
    target += -0.5 * N_g[g] * log_2pi_sigma2 - SS * inv_2sigma2;
  }

  // === Likelihood: Subjects with Missing X ===
  // Marginalize using the fixed max_x from data
  for (k in 1:G_mis) {
    int g = index_mis[k];

    // Compute log probabilities for x = 0 to max_x
    vector[max_x + 1] log_probs;

    for (xg in 0:max_x) {
      real x_val = xg;

      // log p(x_g=k | Z_g)
      real lp_x = neg_binomial_2_lpmf(xg | mu_imp[g], phi);

      // log p(y_g | x_g=k, ...)
      real SS = A_g[g] - 2 * B_g[g] * x_val + C_g[g] * square(x_val);
      real lp_y = -0.5 * N_g[g] * log_2pi_sigma2 - SS * inv_2sigma2;

      log_probs[xg + 1] = lp_x + lp_y;
    }

    target += log_sum_exp(log_probs);
  }
}

generated quantities {
  corr_matrix[2] corr_rand_effects = multiply_lower_tri_self_transpose(L_re);
  cov_matrix[2] cov_rand_effects = quad_form_diag(corr_rand_effects, sigma_re);
}