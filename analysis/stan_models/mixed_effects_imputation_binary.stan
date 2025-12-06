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

  // x is binary and observed only for a subset of subjects
  array[G_obs] int<lower=0, upper=1> x_obs;

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

  // Imputation model (logistic regression for binary x)
  real alpha_imputation;
  vector[S] gamma;
}

transformed parameters {
  // random effects on original scale
  matrix[G, 2] re;
  re = (diag_pre_multiply(sigma_re, L_re) * z_re)';
}

model {
  // Random effects priors
  to_vector(z_re) ~ std_normal();
  sigma_re ~ exponential(0.25);
  L_re ~ lkj_corr_cholesky(2);

  // Main model priors
  sigma_main ~ exponential(0.25);
  alpha_main ~ normal(0, 10);
  beta ~ normal(0, 10);
  beta_t ~ normal(0, 10);
  beta_x ~ normal(0, 10);
  beta_x_t_interaction ~ normal(0, 10);

  // Imputation (logistic) model priors
  alpha_imputation ~ normal(0, 10);
  gamma ~ normal(0, 10);

  // Logistic linear predictor for each subject
  {
    vector[G] eta_imp = alpha_imputation + Z * gamma;

    // ---- Subjects with observed x -----------------------------------------
    for (j in 1:G_obs) {
      int g = index_obs[j];
      int xg = x_obs[j];

      // Logistic regression part for observed x
      target += bernoulli_logit_lpmf(xg | eta_imp[g]);

      // Outcome model for this subject
      for (n in 1:N) {
        if (id[n] == g) {
          real mu = alpha_main
                    + beta_x * xg
                    + re[g, 1]
                    + (beta_t + re[g, 2]) * t[n]
                    + beta_x_t_interaction * (xg * t[n])
                    + X[n] * beta;   // row_vector * vector -> scalar

          target += normal_lpdf(y[n] | mu, sigma_main);
        }
      }
    }

    // ---- Subjects with missing x: marginalize over x in {0, 1} ------------
    for (j in 1:G_mis) {
      int g = index_mis[j];
      real eta_g = eta_imp[g];

      // log p(x_g = 0 | Z_g, params)
      real lp0 = bernoulli_logit_lpmf(0 | eta_g);
      // log p(x_g = 1 | Z_g, params)
      real lp1 = bernoulli_logit_lpmf(1 | eta_g);

      // Add outcome contributions under x=0 and x=1
      for (n in 1:N) {
        if (id[n] == g) {
          // mean if x_g = 0
          real mu0 = alpha_main
                     + 0                         // beta_x * 0
                     + re[g, 1]
                     + (beta_t + re[g, 2]) * t[n]
                     + 0                         // beta_x_t_interaction * 0 * t[n]
                     + X[n] * beta;

          // mean if x_g = 1
          real mu1 = alpha_main
                     + beta_x * 1
                     + re[g, 1]
                     + (beta_t + re[g, 2]) * t[n]
                     + beta_x_t_interaction * t[n]
                     + X[n] * beta;

          lp0 += normal_lpdf(y[n] | mu0, sigma_main);
          lp1 += normal_lpdf(y[n] | mu1, sigma_main);
        }
      }

      // Marginalize over x_g:
      // log p(y_g | Z_g, params) = log [ p(x=0)*p(y_g|x=0) + p(x=1)*p(y_g|x=1) ]
      target += log_sum_exp(lp0, lp1);
    }
  }
}

generated quantities {
  corr_matrix[2] corr_rand_effects =
    multiply_lower_tri_self_transpose(L_re);

  cov_matrix[2] cov_rand_effects =
    quad_form_diag(corr_rand_effects, sigma_re);
}
