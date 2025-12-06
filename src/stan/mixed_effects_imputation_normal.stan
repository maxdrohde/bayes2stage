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

  // Non-centered random effects
  matrix[2, G] z_re;                      // latent standard-normal random effects
  vector<lower=0>[2] sigma_re;            // SDs of random effects
  cholesky_factor_corr[2] L_re;           // Cholesky of random effects correlation

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

  // Missing x values
  vector[G_mis] x_mis;
}

model {
  to_vector(z_re) ~ std_normal();
  sigma_re ~ exponential(0.1);
  L_re ~ lkj_corr_cholesky(2);

  // Main model priors
  sigma_main ~ exponential(0.1);
  alpha_main ~ normal(0, 2);
  beta ~ normal(0, 2);
  beta_t ~ normal(0, 2);
  beta_x ~ normal(0, 2);
  beta_x_t_interaction ~ normal(0, 2);

  // Imputation model priors
  sigma_imputation ~ exponential(0.1);
  alpha_imputation ~ normal(0, 2);
  gamma ~ normal(0, 2);

  vector[G] x;
  x[index_obs] = x_obs;
  x[index_mis] = x_mis;

  matrix[2, G] re;
  re = diag_pre_multiply(sigma_re, L_re) * z_re;

  // Main model
  {
    vector[N] x_subj   = x[id];        // x per observation
    
    vector[N] re_int   = re[1, id]';   // random intercept per obs
    vector[N] re_slope = re[2, id]';   // random slope per obs

    vector[N] alpha_vec =
        alpha_main
      + beta_x * x_subj
      + re_int
      + (beta_t + re_slope) .* t
      + beta_x_t_interaction * (x_subj .* t);

    y ~ normal_id_glm(X, alpha_vec, beta, sigma_main);
  }

  // Imputation model
  {
    x ~ normal_id_glm(Z, alpha_imputation, gamma, sigma_imputation);
  }
}

generated quantities {
  corr_matrix[2] corr_rand_effects =
    multiply_lower_tri_self_transpose(L_re);

  cov_matrix[2] cov_rand_effects =
    quad_form_diag(corr_rand_effects, sigma_re);
}
