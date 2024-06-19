// Assume that there is only one continious covariate with missingness -- denoted x_e
// The other covariates are all at the subject-level (rather than the observation level)

data {
  // Total number of data points
  int<lower=0> N;
  // Number of subjects
  int<lower=0> G;
  // Number of subjects with observed data
  int<lower=0, upper=G> G_obs;
  // Number of subjects with missing data
  int<lower=0, upper=G> G_mis;

  // Number of observed main model predictors
  int<lower=0> P; 

  // Time
  vector[N] t; 
  // Main model covariates
  matrix[N, P] X; 

  // Number of imputation model predictors
  int<lower=0> S;
  // Imputation model covariates
  matrix[G, S] Z;  

  // Outcome
  vector[N] y;

  // Position of observed x_e values
  array[G_obs] int index_obs;
  // Position of missing x_e values
  array[G_mis] int index_mis;
  // Observed x_e values
  vector[G_obs] x_e_obs;

  // Subject ID for each row
  array[N] int<lower=1, upper=G> id; 
}

parameters {
  // Subject-specific random intercepts and slopes
  matrix[G, 2] rand_effects;
  vector<lower=0>[2] sigma_rand_effects;
  corr_matrix[2] corr_rand_effects;
  
  // Main regression parameters
  real alpha_main;
  real beta_t;
  real beta_x_e;
  real beta_xe_t_interaction;
  vector[P] beta;
  real<lower=0> sigma_main;

  // Imputation model parameters
  real alpha_imp;
  vector[S] gamma;
  real<lower=0> sigma_imp;

  // Missing data parameters
  vector[G_mis] x_e_mis;
}

transformed parameters {
  // Fill in the x_e parameter vector with the observed and missing elements
  vector[G] x_e;
  x_e[index_obs] = x_e_obs;
  x_e[index_mis] = x_e_mis;
}

model {
  sigma_main ~ exponential(0.25);          // Prior for main model error SD
  sigma_imp ~ exponential(0.25);           // Prior for imputation model error SD 
  sigma_rand_effects ~ exponential(0.25);  // Prior for random effects SD
  corr_rand_effects ~ lkj_corr(2);         // Prior for random effects correlation

  // Prior for regression parameters
  beta ~ normal(0,10);
  beta_t ~ normal(0,10);
  beta_x_e ~ normal(0,10);
  beta_xe_t_interaction ~ normal(0,10);
  gamma ~ normal(0,10);

  // Model for random intercepts and random slopes
  for (g in 1:G) {
     rand_effects[g,] ~ multi_normal(rep_vector(0,2), quad_form_diag(corr_rand_effects, sigma_rand_effects));
  }
  
  vector[N] intercept;
  for (n in 1:N) {
    intercept[n] =
    alpha_main +
    (beta_x_e * x_e[id[n]]) +
    rand_effects[id[n], 1] +
    (beta_t + rand_effects[id[n], 2]) * t[n] +
    (beta_xe_t_interaction * x_e[id[n]] * t[n]);
  }

  // Main model likelihood
  y ~ normal_id_glm(X, intercept, beta, sigma_main);

  // Imputation model likelihood
  x_e ~ normal_id_glm(Z, alpha_imp, gamma, sigma_imp);
}
