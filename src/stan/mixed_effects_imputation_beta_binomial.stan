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

  // Beta-binomial: number of trials per subject
  array[G] int<lower=1> n_trials;

  // Count covariate for subjects with observed x (values in 0..n_trials[g])
  array[G_obs] int<lower=0> x_obs;

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
  real beta_x;                       // effect of x (raw count)
  real beta_x_t_interaction;
  vector[P] beta;
  real<lower=0> sigma_main;

  // Imputation model for beta-binomial x
  real alpha_imputation;
  vector[S] gamma;
  real<lower=0> phi;                 // dispersion (concentration) parameter
}

model {
  // === Priors ===
  // Random effects
  to_vector(z_re) ~ std_normal();
  sigma_re ~ exponential(0.1);
  L_re ~ lkj_corr_cholesky(2);

  // Main model
  sigma_main ~ exponential(0.1);
  alpha_main ~ normal(0, 100);
  beta ~ normal(0, 100);
  beta_t ~ normal(0, 100);
  beta_x ~ normal(0, 100);
  beta_x_t_interaction ~ normal(0, 100);

  // Imputation model (logit scale)
  alpha_imputation ~ normal(0, 2.5);
  gamma ~ normal(0, 2.5);
  phi ~ exponential(0.1);

  // === Transformed quantities (Local) ===
  
  // Random effects (non-centered)
  matrix[2, G] re = diag_pre_multiply(sigma_re, L_re) * z_re;

  // Mean probability for beta-binomial imputation model
  vector[G] mu_imp = inv_logit(alpha_imputation + Z * gamma);

  // Fixed effects part X*beta
  vector[N] xb = X * beta;

  // === Precompute sufficient statistics for efficient marginalization ===
  // We expand the squared error: (resid - x_val * delta)^2
  // into: resid^2 - 2*x_val*(resid*delta) + x_val^2*(delta^2)
  
  vector[G] N_g = rep_vector(0, G);   // observation count per subject
  vector[G] A_g = rep_vector(0, G);   // sum of squared base residuals
  vector[G] B_g = rep_vector(0, G);   // sum of residual * delta
  vector[G] C_g = rep_vector(0, G);   // sum of delta^2

  for (n in 1:N) {
    int g = id[n];
    
    // Base expectation (assuming x=0)
    real base = alpha_main + re[1, g] + (beta_t + re[2, g]) * t[n] + xb[n];
    real resid = y[n] - base;
    
    // Slope contribution of x
    real delta = beta_x + beta_x_t_interaction * t[n];

    N_g[g] += 1;
    A_g[g] += square(resid);
    B_g[g] += resid * delta;
    C_g[g] += square(delta);
  }

  // Precompute likelihood constants
  real log_2pi_sigma2 = log(2 * pi() * square(sigma_main));
  real inv_2sigma2 = 0.5 / square(sigma_main);

  // === Likelihood: Subjects with Observed X ===
  for (k in 1:G_obs) {
    int g = index_obs[k];
    int n_g = n_trials[g];
    int xg = x_obs[k];
    
    // UPDATED: Use raw count xg, not proportion
    real x_val = xg; 

    // Beta-binomial parameters (alpha/beta count formulation)
    real alpha_bb = mu_imp[g] * phi;
    real beta_bb = (1 - mu_imp[g]) * phi;

    // 1. Imputation model: p(x_g | Z_g)
    target += beta_binomial_lpmf(xg | n_g, alpha_bb, beta_bb);

    // 2. Main model: p(y_g | x_g, ...) using sufficient stats
    real SS = A_g[g] - 2 * B_g[g] * x_val + C_g[g] * square(x_val);
    target += -0.5 * N_g[g] * log_2pi_sigma2 - SS * inv_2sigma2;
  }

  // === Likelihood: Subjects with Missing X ===
  // Marginalize over discrete x_g in {0, 1, ..., n_trials_g}
  for (k in 1:G_mis) {
    int g = index_mis[k];
    int n_g = n_trials[g];

    real alpha_bb = mu_imp[g] * phi;
    real beta_bb = (1 - mu_imp[g]) * phi;

    // Vector to hold log probability for each possible integer value of x
    vector[n_g + 1] log_probs;

    for (xg in 0:n_g) {
      // UPDATED: Use raw count xg, not proportion
      real x_val = xg;

      // log p(x_g=k | Z_g)
      real lp_x = beta_binomial_lpmf(xg | n_g, alpha_bb, beta_bb);

      // log p(y_g | x_g=k, ...)
      real SS = A_g[g] - 2 * B_g[g] * x_val + C_g[g] * square(x_val);
      real lp_y = -0.5 * N_g[g] * log_2pi_sigma2 - SS * inv_2sigma2;

      log_probs[xg + 1] = lp_x + lp_y;
    }

    // Marginalize: log sum( exp(log_probs) )
    target += log_sum_exp(log_probs);
  }
}

generated quantities {
  // Random Effects Correlation and Covariance matrices
  corr_matrix[2] corr_rand_effects = multiply_lower_tri_self_transpose(L_re);
  cov_matrix[2] cov_rand_effects = quad_form_diag(corr_rand_effects, sigma_re);
}





