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

  // binary subject-level covariate for subjects with observed x
  array[G_obs] int<lower=0, upper=1> x_obs;

  // subject index per row
  array[N] int<lower=1, upper=G> id;
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

  // Imputation model for binary x: Bernoulli-logit(alpha_imputation + Z * gamma)
  real alpha_imputation;
  vector[S] gamma;
}

model {
  // Random effects priors
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

  // Imputation model priors (logit scale, following Gelman et al. recommendations)
  alpha_imputation ~ normal(0, 2.5);
  gamma ~ normal(0, 2.5);

  // Random effects (non-centered parameterization)
  matrix[2, G] re = diag_pre_multiply(sigma_re, L_re) * z_re;

  // Imputation model linear predictor for subject-level binary x
  vector[G] eta_imputation = alpha_imputation + Z * gamma;

  // Precompute per-subject log-likelihood of y given x_g = 0 and x_g = 1
  vector[G] log_y_x0 = rep_vector(0, G);
  vector[G] log_y_x1 = rep_vector(0, G);

  // Contribution of X * beta reused across both x = 0 and x = 1 cases
  vector[N] xb = X * beta;

  // Build subject-wise log-likelihood pieces
  for (n in 1:N) {
    int g = id[n];

    // Part of the linear predictor that does *not* depend on x_g
    real base =
      alpha_main
      + re[1, g]                               // random intercept
      + (beta_t + re[2, g]) * t[n]            // time slope + random slope
      + xb[n];                                // fixed effects X*beta

    // Extra contribution when x_g = 1 (vs x_g = 0)
    real delta = beta_x + beta_x_t_interaction * t[n];

    // Log-likelihood contributions for this observation under x_g = 0 and x_g = 1
    log_y_x0[g] += normal_lpdf(y[n] | base,          sigma_main);
    log_y_x1[g] += normal_lpdf(y[n] | base + delta,  sigma_main);
  }

  // Subjects with observed binary x: explicit Bernoulli-logit likelihood + y|x
  for (k in 1:G_obs) {
    int g  = index_obs[k];
    int xg = x_obs[k];

    // Imputation model: p(x_g | Z_g)
    target += bernoulli_logit_lpmf(xg | eta_imputation[g]);

    // Main model: p(y_g | x_g, ...)
    if (xg == 1) {
      target += log_y_x1[g];
    } else {
      target += log_y_x0[g];
    }
  }

  // Subjects with missing x: marginalize x_g ∈ {0, 1}
  //
  // log p(y_g | Z_g, ...) =
  //   log_sum_exp( log p(x_g=1|Z_g) + log p(y_g|x_g=1),
  //                log p(x_g=0|Z_g) + log p(y_g|x_g=0) )
  //
  // Using bernoulli_logit_lpmf to stay on the logit scale for numerical stability.
  for (k in 1:G_mis) {
    int g = index_mis[k];

    real lp_x1 = bernoulli_logit_lpmf(1 | eta_imputation[g]);
    real lp_x0 = bernoulli_logit_lpmf(0 | eta_imputation[g]);

    real lp1 = lp_x1 + log_y_x1[g];
    real lp0 = lp_x0 + log_y_x0[g];

    target += log_sum_exp(lp1, lp0);
  }
}

generated quantities {
  corr_matrix[2] corr_rand_effects =
    multiply_lower_tri_self_transpose(L_re);

  cov_matrix[2] cov_rand_effects =
    quad_form_diag(corr_rand_effects, sigma_re);
}
