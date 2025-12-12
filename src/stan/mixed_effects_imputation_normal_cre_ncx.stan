// Hybrid: Centered random effects + Non-centered x imputation
// Best for: Large N (well-identified RE) with sparse x sampling
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
  vector[G_obs] x_obs;

  array[N] int<lower=1, upper=G> id;  // subject index per row
}

parameters {
  // Centered random effects (good for large N)
  matrix[G, 2] re;                      // random intercepts and slopes
  vector<lower=0>[2] sigma_re;          // SDs of random effects
  cholesky_factor_corr[2] L_re;         // Cholesky of random effects correlation

  // Main model
  real alpha_main;
  real beta_t;
  real beta_x;
  real beta_x_t_interaction;
  vector[P] beta;
  real<lower=0> sigma_main;

  // Imputation model
  real alpha_imputation;
  vector[S] gamma;
  real<lower=0> sigma_imputation;

  // Non-centered missing x values
  vector[G_mis] x_mis_raw;
}

transformed parameters {
  // Imputation model linear predictor
  vector[G] x_imputation_lp = alpha_imputation + Z * gamma;

  // Non-centered parameterization for missing x
  vector[G_mis] x_mis = x_imputation_lp[index_mis] + sigma_imputation * x_mis_raw;

  // Full x vector (observed + imputed)
  vector[G] x;
  x[index_obs] = x_obs;
  x[index_mis] = x_mis;
}

model {
  // Random effects hyperpriors
  sigma_re ~ exponential(0.1);
  L_re ~ lkj_corr_cholesky(2);

  // Centered random effects prior
  {
    matrix[2, 2] L_Sigma = diag_pre_multiply(sigma_re, L_re);
    for (g in 1:G) {
      re[g]' ~ multi_normal_cholesky(rep_vector(0, 2), L_Sigma);
    }
  }

  // Main model priors
  sigma_main ~ exponential(0.1);
  alpha_main ~ normal(0, 100);
  beta ~ normal(0, 100);
  beta_t ~ normal(0, 100);
  beta_x ~ normal(0, 100);
  beta_x_t_interaction ~ normal(0, 100);

  // Imputation model priors
  sigma_imputation ~ exponential(0.1);
  alpha_imputation ~ normal(0, 100);
  gamma ~ normal(0, 100);

  // Non-centered prior for missing x
  x_mis_raw ~ std_normal();

  // Imputation model likelihood (observed x only)
  x_obs ~ normal(x_imputation_lp[index_obs], sigma_imputation);

  // Main model likelihood
  {
    vector[N] x_subj   = x[id];        // x per observation
    vector[N] re_int   = re[id, 1];    // random intercept per obs
    vector[N] re_slope = re[id, 2];    // random slope per obs

    vector[N] alpha_vec =
        alpha_main
      + beta_x * x_subj
      + re_int
      + (beta_t + re_slope) .* t
      + beta_x_t_interaction * (x_subj .* t);

    y ~ normal_id_glm(X, alpha_vec, beta, sigma_main);
  }
}

generated quantities {
  corr_matrix[2] corr_rand_effects =
    multiply_lower_tri_self_transpose(L_re);

  cov_matrix[2] cov_rand_effects =
    quad_form_diag(corr_rand_effects, sigma_re);
}
