// Mixed Effects Model with Marginalized Random Effects AND Marginalized Binary X
//
// This model integrates out BOTH:
//   1. Random effects (analytically via MVN)
//   2. Missing binary x values (via mixture of two Gaussians)
//
// Result: ZERO latent variables for subjects with missing x
//
// For subjects with missing x:
//   p(y_g | theta) = p(y_g | x=0, theta) * p(x=0 | theta) + p(y_g | x=1, theta) * p(x=1 | theta)
//
// where p(y_g | x, theta) is MVN after marginalizing out random effects.
//
// Trade-off:
//   - More expensive likelihood (2 MVN densities per subject with missing x)
//   - But eliminates ALL latent variables (no funnel, no discrete sampling issues)
//
// Best for: Large G with binary x causing mixing problems

functions {
    /**
     * Compute marginal covariance matrix for a subject
     *
     * V_g = Z_re * Sigma_re * Z_re' + sigma^2 * I
     *
     * where Z_re = [1, t_g] is n_g x 2
     *
     * @param t_g Vector of time points (length n_g)
     * @param L_Sigma Lower Cholesky factor of RE covariance (2x2)
     * @param sigma Residual standard deviation
     * @return Cholesky factor of marginal covariance (n_g x n_g)
     */
    matrix compute_marginal_cov_chol(vector t_g, matrix L_Sigma, real sigma) {
        int n = rows(t_g);
        matrix[n, 2] Z_re;
        matrix[n, n] V;

        // Build RE design matrix: [1, t]
        Z_re[, 1] = rep_vector(1.0, n);
        Z_re[, 2] = t_g;

        // V = Z * L * L' * Z' + sigma^2 * I
        V = tcrossprod(Z_re * L_Sigma);

        // Add residual variance to diagonal
        for (i in 1:n) {
            V[i, i] += square(sigma);
        }

        return cholesky_decompose(V);
    }

    /**
     * Compute mean vector for subject g given x value
     *
     * @param x Value of binary covariate (0 or 1)
     * @param t_g Time vector for subject
     * @param X_g Covariate matrix for subject (n_g x P)
     * @param alpha_main Intercept
     * @param beta_x Effect of x
     * @param beta_t Effect of time
     * @param beta_x_t x-time interaction
     * @param beta Other covariate effects
     * @return Mean vector (n_g)
     */
    vector compute_mean(int x, vector t_g, matrix X_g,
                        real alpha_main, real beta_x, real beta_t,
                        real beta_x_t, vector beta) {
        int n = rows(t_g);
        vector[n] mu = rep_vector(alpha_main + beta_x * x, n)
                       + beta_t * t_g
                       + beta_x_t * x * t_g
                       + X_g * beta;
        return mu;
    }
}

data {
    // Dimensions
    int<lower=0> N;                         // Total observations
    int<lower=1> G;                         // Total subjects
    int<lower=0, upper=G> G_obs;            // Subjects with observed x
    int<lower=0, upper=G> G_mis;            // Subjects with missing x

    int<lower=0> P;                         // Main model predictors (excluding x)
    vector[N] t;                            // Time variable
    matrix[N, P] X;                         // Main model covariates

    int<lower=0> S;                         // Imputation model predictors
    matrix[G, S] Z;                         // Imputation covariates (subject-level)

    vector[N] y;                            // Outcome variable

    // Subject indexing
    array[G_obs] int<lower=1, upper=G> index_obs;  // Subjects with observed x
    array[G_mis] int<lower=1, upper=G> index_mis;  // Subjects with missing x
    array[G_obs] int<lower=0, upper=1> x_obs;      // Observed x values (binary)

    array[N] int<lower=1, upper=G> id;      // Subject index per observation

    // Subject start position and length (for ragged array access)
    array[G] int<lower=1> pos;              // Start index for subject g
    array[G] int<lower=1> len;              // Number of obs for subject g
}

parameters {
    // Random effects variance components (marginalized, but still estimated)
    vector<lower=0>[2] sigma_re;            // [intercept SD, slope SD]
    cholesky_factor_corr[2] L_re;           // Cholesky factor of RE correlation

    // Main model fixed effects
    real alpha_main;
    real beta_t;
    real beta_x;
    real beta_x_t_interaction;
    vector[P] beta;
    real<lower=0> sigma_main;

    // Imputation model (Bernoulli logistic regression)
    real alpha_imputation;
    vector[S] gamma;

    // NOTE: No x_mis parameters! Binary x is marginalized out.
    // NOTE: No z_re parameters! Random effects are marginalized out.
}

transformed parameters {
    // Cholesky factor of RE covariance: L_Sigma = diag(sigma_re) * L_corr
    matrix[2, 2] L_Sigma = diag_pre_multiply(sigma_re, L_re);

    // Imputation model linear predictor for all subjects
    vector[G] eta_imputation = alpha_imputation + Z * gamma;
}

model {
    // =========================================================================
    // Priors
    // =========================================================================

    // Random effects variance components
    sigma_re ~ exponential(0.1);
    L_re ~ lkj_corr_cholesky(2);

    // Main model priors
    sigma_main ~ exponential(0.1);
    alpha_main ~ normal(0, 100);
    beta ~ normal(0, 100);
    beta_t ~ normal(0, 100);
    beta_x ~ normal(0, 100);
    beta_x_t_interaction ~ normal(0, 100);

    // Imputation model priors
    alpha_imputation ~ normal(0, 2.5);
    gamma ~ normal(0, 2.5);

    // =========================================================================
    // Likelihood for Subjects with OBSERVED x
    // =========================================================================
    //
    // For these subjects:
    //   p(y_g, x_g | theta) = p(y_g | x_g, theta) * p(x_g | theta)
    //
    // where p(y_g | x_g, theta) is MVN after marginalizing RE

    for (k in 1:G_obs) {
        int g = index_obs[k];
        int x_g = x_obs[k];
        int n_g = len[g];
        int start_idx = pos[g];

        // Extract subject's data
        vector[n_g] y_g = segment(y, start_idx, n_g);
        vector[n_g] t_g = segment(t, start_idx, n_g);
        matrix[n_g, P] X_g = X[start_idx:(start_idx + n_g - 1), ];

        // Compute marginal covariance (Cholesky factor)
        matrix[n_g, n_g] L_V = compute_marginal_cov_chol(t_g, L_Sigma, sigma_main);

        // Compute mean given observed x
        vector[n_g] mu_g = compute_mean(x_g, t_g, X_g,
                                        alpha_main, beta_x, beta_t,
                                        beta_x_t_interaction, beta);

        // Outcome likelihood (marginalized over RE)
        y_g ~ multi_normal_cholesky(mu_g, L_V);

        // Imputation model likelihood
        x_g ~ bernoulli_logit(eta_imputation[g]);
    }

    // =========================================================================
    // Likelihood for Subjects with MISSING x (Fully Marginalized)
    // =========================================================================
    //
    // For these subjects:
    //   p(y_g | theta) = p(y_g | x=0, theta) * p(x=0 | theta)
    //                  + p(y_g | x=1, theta) * p(x=1 | theta)
    //
    // This is a mixture of two MVN distributions.
    // We use log_sum_exp for numerical stability.

    for (k in 1:G_mis) {
        int g = index_mis[k];
        int n_g = len[g];
        int start_idx = pos[g];

        // Extract subject's data
        vector[n_g] y_g = segment(y, start_idx, n_g);
        vector[n_g] t_g = segment(t, start_idx, n_g);
        matrix[n_g, P] X_g = X[start_idx:(start_idx + n_g - 1), ];

        // Compute marginal covariance (same for x=0 and x=1)
        matrix[n_g, n_g] L_V = compute_marginal_cov_chol(t_g, L_Sigma, sigma_main);

        // Compute mean for x = 0
        vector[n_g] mu_0 = compute_mean(0, t_g, X_g,
                                        alpha_main, beta_x, beta_t,
                                        beta_x_t_interaction, beta);

        // Compute mean for x = 1
        vector[n_g] mu_1 = compute_mean(1, t_g, X_g,
                                        alpha_main, beta_x, beta_t,
                                        beta_x_t_interaction, beta);

        // Log densities for outcome given x
        real lp_y_x0 = multi_normal_cholesky_lpdf(y_g | mu_0, L_V);
        real lp_y_x1 = multi_normal_cholesky_lpdf(y_g | mu_1, L_V);

        // Log probabilities for x from imputation model
        // log(p) = log(inv_logit(eta)) = -log1p_exp(-eta)
        // log(1-p) = log(1 - inv_logit(eta)) = -log1p_exp(eta)
        real lp_x1 = -log1p_exp(-eta_imputation[g]);  // log p(x=1)
        real lp_x0 = -log1p_exp(eta_imputation[g]);   // log p(x=0)

        // Marginal likelihood via log-sum-exp
        target += log_sum_exp(lp_y_x0 + lp_x0, lp_y_x1 + lp_x1);
    }
}

generated quantities {
    // Correlation matrix for random effects
    corr_matrix[2] corr_rand_effects = multiply_lower_tri_self_transpose(L_re);

    // Covariance matrix for random effects
    cov_matrix[2] cov_rand_effects = quad_form_diag(corr_rand_effects, sigma_re);

    // Log likelihood for LOO-CV (per subject)
    vector[G] log_lik;

    // Compute log likelihoods for observed subjects
    for (k in 1:G_obs) {
        int g = index_obs[k];
        int x_g = x_obs[k];
        int n_g = len[g];
        int start_idx = pos[g];

        vector[n_g] y_g = segment(y, start_idx, n_g);
        vector[n_g] t_g = segment(t, start_idx, n_g);
        matrix[n_g, P] X_g = X[start_idx:(start_idx + n_g - 1), ];

        matrix[n_g, n_g] L_V = compute_marginal_cov_chol(t_g, L_Sigma, sigma_main);
        vector[n_g] mu_g = compute_mean(x_g, t_g, X_g,
                                        alpha_main, beta_x, beta_t,
                                        beta_x_t_interaction, beta);

        real lp_y = multi_normal_cholesky_lpdf(y_g | mu_g, L_V);
        real lp_x = bernoulli_logit_lpmf(x_g | eta_imputation[g]);

        log_lik[g] = lp_y + lp_x;
    }

    // Compute log likelihoods for missing subjects
    for (k in 1:G_mis) {
        int g = index_mis[k];
        int n_g = len[g];
        int start_idx = pos[g];

        vector[n_g] y_g = segment(y, start_idx, n_g);
        vector[n_g] t_g = segment(t, start_idx, n_g);
        matrix[n_g, P] X_g = X[start_idx:(start_idx + n_g - 1), ];

        matrix[n_g, n_g] L_V = compute_marginal_cov_chol(t_g, L_Sigma, sigma_main);

        vector[n_g] mu_0 = compute_mean(0, t_g, X_g,
                                        alpha_main, beta_x, beta_t,
                                        beta_x_t_interaction, beta);
        vector[n_g] mu_1 = compute_mean(1, t_g, X_g,
                                        alpha_main, beta_x, beta_t,
                                        beta_x_t_interaction, beta);

        real lp_y_x0 = multi_normal_cholesky_lpdf(y_g | mu_0, L_V);
        real lp_y_x1 = multi_normal_cholesky_lpdf(y_g | mu_1, L_V);

        real lp_x1 = -log1p_exp(-eta_imputation[g]);
        real lp_x0 = -log1p_exp(eta_imputation[g]);

        // Marginal log likelihood
        log_lik[g] = log_sum_exp(lp_y_x0 + lp_x0, lp_y_x1 + lp_x1);
    }
}
