data {
  int<lower=0> N;
  int<lower=1> G;
  int<lower=0, upper=G> G_obs;
  int<lower=0, upper=G> G_mis;

  int<lower=0> P;
  vector[N] t;
  matrix[N, P] X;

  int<lower=0> S;
  matrix[G, S] Z;

  vector[N] y;

  array[G_obs] int<lower=1, upper=G> index_obs;
  array[G_mis] int<lower=1, upper=G> index_mis;

  array[G_obs] int<lower=0> x_obs;
  
  int<lower=max(x_obs)> max_x;

  // IMPORTANT: Data must be sorted by ID for this optimization to work efficiently
  array[N] int<lower=1, upper=G> id;
}

transformed data {
  // Precompute subject indices to allow vectorization
  // This assumes the data is sorted by ID. 
  array[G] int pos_start;
  array[G] int pos_end;
  array[G] int len;
  
  {
    int current_g = 0;
    for (n in 1:N) {
      if (id[n] != current_g) {
        current_g = id[n];
        pos_start[current_g] = n;
      }
      pos_end[current_g] = n;
    }
    for (g in 1:G) {
      len[g] = pos_end[g] - pos_start[g] + 1;
    }
  }

  // Precompute sequence for vectorization
  vector[max_x + 1] x_seq;
  vector[max_x + 1] x_seq_sq;
  for (i in 0:max_x) {
    x_seq[i+1] = i;
    x_seq_sq[i+1] = square(i);
  }
}

parameters {
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
  real<lower=0> phi;                 
}

model {
  // === Priors ===
  to_vector(z_re) ~ std_normal();
  sigma_re ~ exponential(0.1);
  L_re ~ lkj_corr_cholesky(2);

  sigma_main ~ exponential(0.1);
  alpha_main ~ normal(0, 100);
  beta ~ normal(0, 100);
  beta_t ~ normal(0, 100);
  beta_x ~ normal(0, 100);
  beta_x_t_interaction ~ normal(0, 100);

  // Imputation model (log scale)
  alpha_imputation ~ normal(0, 2.5);
  gamma ~ normal(0, 2.5);
  phi ~ exponential(0.1);

  // === Transformed quantities (Local) ===
  matrix[2, G] re = diag_pre_multiply(sigma_re, L_re) * z_re;
  vector[G] mu_imp = exp(alpha_imputation + Z * gamma);
  vector[N] xb = X * beta; // Global matrix multiplication is efficient

  // Constants for likelihood
  real log_2pi_sigma2 = log(2 * pi() * square(sigma_main));
  real inv_2sigma2 = 0.5 / square(sigma_main);

  // === Main Calculation Loop (Subject-wise) ===
  // We compute sufficient stats and likelihood in one pass over G
  
  // Buffers for sufficient stats
  real Ag; 
  real Bg; 
  real Cg;
  
  // Loop over all subjects to compute sufficient stats efficiently
  // We store them temporarily to use in the Obs/Mis blocks below
  // Note: To further optimize, you could split this into two loops (one for Obs, one for Mis)
  // to avoid storing Ag, Bg, Cg in vectors, but this is already much faster.
  
  vector[G] A_vec;
  vector[G] B_vec;
  vector[G] C_vec;

  for (g in 1:G) {
    // Extract subject data segments
    int s = pos_start[g];
    int l = len[g];
    
    // Vector slices
    vector[l] t_sub = segment(t, s, l);
    vector[l] y_sub = segment(y, s, l);
    vector[l] xb_sub = segment(xb, s, l);

    // Vectorized Base calculation
    // base = alpha + b0 + (beta_t + b1)*t + xb
    vector[l] base = alpha_main + re[1, g] + (beta_t + re[2, g]) * t_sub + xb_sub;
    vector[l] resid = y_sub - base;
    vector[l] delta = beta_x + beta_x_t_interaction * t_sub;

    // Vectorized Sufficient Statistics (Dot products are fast)
    A_vec[g] = dot_self(resid);        // Sum of squares
    B_vec[g] = dot_product(resid, delta);
    C_vec[g] = dot_self(delta);
  }

  // === Likelihood: Subjects with Observed X ===
  for (k in 1:G_obs) {
    int g = index_obs[k];
    int xg = x_obs[k];
    
    // 1. Imputation model
    target += neg_binomial_2_lpmf(xg | mu_imp[g], phi);

    // 2. Main model
    real SS = A_vec[g] - 2 * B_vec[g] * xg + C_vec[g] * square(xg);
    target += -0.5 * len[g] * log_2pi_sigma2 - SS * inv_2sigma2;
  }

  // === Likelihood: Subjects with Missing X ===
  // OPTIMIZATION: Recurrence relation + Vectorization
  for (k in 1:G_mis) {
    int g = index_mis[k];
    real mu = mu_imp[g];
    
    // 1. Calculate Log-Probabilities for NB via Recurrence
    // log P(0) = phi * log(phi / (mu + phi))
    // P(k) = P(k-1) * (k - 1 + phi) / k * (mu / (mu + phi))
    
    vector[max_x + 1] log_probs_nb;
    real log_r_mu = log(mu / (mu + phi));
    
    // Base case (x=0)
    log_probs_nb[1] = phi * log(phi / (mu + phi)); 

    // Recurrence loop (x=1 to max_x)
    for (i in 1:max_x) {
       log_probs_nb[i+1] = log_probs_nb[i] + log(i - 1.0 + phi) - log(i) + log_r_mu;
    }

    // 2. Calculate Log-Likelihood for Y (Vectorized)
    // SS = A - 2Bx + Cx^2
    vector[max_x + 1] SS_vec = A_vec[g] - 2 * B_vec[g] * x_seq + C_vec[g] * x_seq_sq;
    
    vector[max_x + 1] lp_y = -0.5 * len[g] * log_2pi_sigma2 - SS_vec * inv_2sigma2;

    // 3. Combine and marginalize
    target += log_sum_exp(log_probs_nb + lp_y);
  }
}

generated quantities {
  corr_matrix[2] corr_rand_effects = multiply_lower_tri_self_transpose(L_re);
  cov_matrix[2] cov_rand_effects = quad_form_diag(corr_rand_effects, sigma_re);
}
