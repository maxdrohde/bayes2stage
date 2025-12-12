// Marginal Mixed Effects Model with Imputation
// =============================================================================
// This model integrates out random effects analytically, dramatically reducing
// the number of parameters from O(G) to O(1) for the random effects component.
//
// For N=2000, this reduces parameters from ~5500 to ~1500 (just missing x values)
// which greatly improves MCMC mixing.
//
// Key features:
// - Random effects are marginalized (integrated out analytically)
// - Uses Cholesky decomposition for efficient covariance computation
// - Weakly informative priors with configurable scales
// - Handles missing x via joint imputation model
// =============================================================================

data {
  int<lower=0> N;                   // total number of observations
  int<lower=1> G;                   // number of subjects
  int<lower=0, upper=G> G_obs;      // number of subjects with observed x
  int<lower=0, upper=G> G_mis;      // number of subjects with missing x

  int<lower=0> P;                   // number of main model covariates
  vector[N] t;                      // time (flat vector)
  matrix[N, P] X;                   // main model covariates

  int<lower=0> S;                   // number of imputation model covariates
  matrix[G, S] Z;                   // imputation model covariates

  vector[N] y;                      // outcome (flat vector)

  array[G_obs] int<lower=1, upper=G> index_obs;
  array[G_mis] int<lower=1, upper=G> index_mis;
  vector[G_obs] x_obs;

  array[N] int<lower=1, upper=G> id;  // subject index per observation

  // Subject observation indices (for grouping)
  array[G] int<lower=1> pos;         // start position for each subject
  array[G] int<lower=1> len;         // number of observations per subject

  // Prior scales (allows tuning for different datasets)
  real<lower=0> prior_beta_sd;          // SD for regression coefficient priors (default: 5)
  real<lower=0> prior_sigma_rate;       // Rate for exponential prior on sigmas (default: 1)
  real<lower=0> prior_sigma_re_rate;    // Rate for exponential prior on RE sigmas (default: 0.5)
}

parameters {
  // Main model fixed effects
  real alpha_main;
  real beta_t;
  real beta_x;
  real beta_x_t_interaction;
  vector[P] beta;
  real<lower=0> sigma_main;

  // Random effects variance components
  vector<lower=0>[2] sigma_re;          // SDs: [intercept, slope]
  cholesky_factor_corr[2] L_corr_re;    // Cholesky factor of correlation matrix

  // Imputation model
  real alpha_imputation;
  vector[S] gamma;
  real<lower=0> sigma_imputation;

  // Missing x values (directly parameterized - centered)
  vector[G_mis] x_mis;
}

transformed parameters {
  // Full x vector (observed + imputed)
  vector[G] x_full;
  x_full[index_obs] = x_obs;
  x_full[index_mis] = x_mis;

  // Covariance matrix for random effects (2x2)
  matrix[2, 2] Sigma_re = diag_pre_multiply(sigma_re, L_corr_re)
                        * diag_pre_multiply(sigma_re, L_corr_re)';
}

model {
  // =========================================================================
  // Priors (weakly informative, appropriately scaled)
  // =========================================================================

  // Fixed effects - normal priors centered at 0
  alpha_main ~ normal(0, prior_beta_sd);
  beta_t ~ normal(0, prior_beta_sd);
  beta_x ~ normal(0, prior_beta_sd);
  beta_x_t_interaction ~ normal(0, prior_beta_sd);
  beta ~ normal(0, prior_beta_sd);

  // Residual SD
  sigma_main ~ exponential(prior_sigma_rate);

  // Random effect SDs
  sigma_re ~ exponential(prior_sigma_re_rate);

  // Correlation matrix - LKJ prior (eta=2 favors moderate correlations)
  L_corr_re ~ lkj_corr_cholesky(2);

  // Imputation model priors
  alpha_imputation ~ normal(0, prior_beta_sd);
  gamma ~ normal(0, prior_beta_sd);
  sigma_imputation ~ exponential(prior_sigma_rate);

  // =========================================================================
  // Imputation model likelihood
  // x_i ~ Normal(alpha_imp + Z_i * gamma, sigma_imp)
  // =========================================================================
  {
    vector[G] mu_x = alpha_imputation + Z * gamma;
    x_full ~ normal(mu_x, sigma_imputation);
  }

  // =========================================================================
  // Marginal likelihood for y (random effects integrated out analytically)
  // =========================================================================
  //
  // Model: y_ij = alpha + beta_x*x_i + X_ij*beta + (beta_t + u_1i)*t_ij + u_0i + eps_ij
  //
  // For subject i with M_i observations:
  //   y_i | x_i ~ Normal(mu_i, V_i)
  //
  // where:
  //   mu_i = alpha + beta_x*x_i + X_i*beta + beta_t*t_i + beta_xt*(x_i*t_i)
  //   V_i = Z_i * Sigma_re * Z_i' + sigma^2 * I
  //   Z_i = [1, t_i]  (M_i x 2 random effects design matrix)
  //
  // This integration reduces ~4000 random effect parameters to just 3 variance
  // components (sigma_re[1], sigma_re[2], and correlation).
  // =========================================================================

  real sigma2_main = square(sigma_main);

  for (g in 1:G) {
    int M_g = len[g];
    int start_g = pos[g];
    int end_g = start_g + M_g - 1;

    // Extract data for subject g
    vector[M_g] y_g = y[start_g:end_g];
    vector[M_g] t_g = t[start_g:end_g];
    matrix[M_g, P] X_g = X[start_g:end_g, ];

    // Random effects design matrix: Z_g = [1, t]
    matrix[M_g, 2] Z_re_g;
    Z_re_g[, 1] = rep_vector(1.0, M_g);
    Z_re_g[, 2] = t_g;

    // Marginal mean
    vector[M_g] mu_g = rep_vector(alpha_main + beta_x * x_full[g], M_g)
                     + X_g * beta
                     + beta_t * t_g
                     + beta_x_t_interaction * x_full[g] * t_g;

    // Marginal covariance: V_g = Z_g * Sigma_re * Z_g' + sigma^2 * I
    matrix[M_g, M_g] V_g = Z_re_g * Sigma_re * Z_re_g'
                        + diag_matrix(rep_vector(sigma2_main, M_g));

    // Add to log probability using Cholesky for numerical stability
    y_g ~ multi_normal(mu_g, V_g);
  }
}

generated quantities {
  // Correlation matrix for random effects
  corr_matrix[2] corr_rand_effects = multiply_lower_tri_self_transpose(L_corr_re);

  // Covariance matrix
  cov_matrix[2] cov_rand_effects = Sigma_re;
}
