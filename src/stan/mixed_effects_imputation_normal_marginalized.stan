// Mixed Effects Model with Fully Marginalized Random Effects AND Missing X
//
// This model integrates out BOTH:
//   1. Random effects (analytically via MVN)
//   2. Missing continuous x values (analytically via Gaussian conjugacy)
//
// Result: ZERO latent variables
//
// For subjects with missing x:
//   If μ(x) = c + x*d and x ~ N(μ_x, σ²_x), then:
//   y ~ MVN(c + μ_x*d, V + σ²_x*d*d')
//
// This is a rank-1 update to the covariance matrix.
//
// Trade-off:
//   - More expensive likelihood (MVN with modified covariance per subject)
//   - But eliminates ALL latent variables (no funnel, perfect mixing)
//
// Best for: Large G with continuous x causing mixing problems

functions {
    /**
     * Compute Cholesky factor of marginal covariance matrix for a subject
     *
     * V_g = Z_re * Sigma_re * Z_re' + sigma^2 * I
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
     * Compute Cholesky factor of marginal covariance with x uncertainty
     *
     * V_aug = V + sigma_x^2 * d * d'
     *
     * where d[i] = beta_x + beta_x_t * t[i] (the coefficient on x)
     *
     * This is a rank-1 update: if V = L*L', then
     * V + sigma_x^2 * d*d' can be computed via Cholesky update
     *
     * @param L_V Cholesky factor of base covariance (n x n)
     * @param d Vector of x coefficients (length n)
     * @param sigma_x Standard deviation of x
     * @return Cholesky factor of augmented covariance (n x n)
     */
    matrix chol_rank1_update(matrix L_V, vector d, real sigma_x) {
        int n = rows(L_V);
        matrix[n, n] V = tcrossprod(L_V);

        // Rank-1 update: V_aug = V + sigma_x^2 * d * d'
        V += square(sigma_x) * tcrossprod(to_matrix(d, n, 1));

        return cholesky_decompose(V);
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
    vector[G_obs] x_obs;                           // Observed x values

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

    // Imputation model
    real alpha_imputation;
    vector[S] gamma;
    real<lower=0> sigma_imputation;

    // NOTE: No x_mis parameters! Continuous x is marginalized out.
    // NOTE: No z_re parameters! Random effects are marginalized out.
}

transformed parameters {
    // Cholesky factor of RE covariance: L_Sigma = diag(sigma_re) * L_corr
    matrix[2, 2] L_Sigma = diag_pre_multiply(sigma_re, L_re);

    // Imputation model linear predictor (E[x] for each subject)
    vector[G] x_imputation_mean = alpha_imputation + Z * gamma;
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
    sigma_imputation ~ exponential(0.1);
    alpha_imputation ~ normal(0, 100);
    gamma ~ normal(0, 100);

    // =========================================================================
    // Likelihood for Subjects with OBSERVED x
    // =========================================================================

    for (k in 1:G_obs) {
        int g = index_obs[k];
        real x_g = x_obs[k];
        int n_g = len[g];
        int start_idx = pos[g];

        // Extract subject's data
        vector[n_g] y_g = segment(y, start_idx, n_g);
        vector[n_g] t_g = segment(t, start_idx, n_g);

        // Compute base mean (without x terms)
        vector[n_g] c_g;
        for (i in 1:n_g) {
            int obs_idx = start_idx + i - 1;
            c_g[i] = alpha_main + beta_t * t_g[i] + dot_product(X[obs_idx], beta);
        }

        // x coefficient vector: d[i] = beta_x + beta_x_t * t[i]
        vector[n_g] d_g = rep_vector(beta_x, n_g) + beta_x_t_interaction * t_g;

        // Full mean: mu = c + x * d
        vector[n_g] mu_g = c_g + x_g * d_g;

        // Marginal covariance (integrating out RE only)
        matrix[n_g, n_g] L_V = compute_marginal_cov_chol(t_g, L_Sigma, sigma_main);

        // Outcome likelihood
        y_g ~ multi_normal_cholesky(mu_g, L_V);

        // Imputation model likelihood for observed x
        x_g ~ normal(x_imputation_mean[g], sigma_imputation);
    }

    // =========================================================================
    // Likelihood for Subjects with MISSING x (Fully Marginalized)
    // =========================================================================
    //
    // For these subjects, we integrate out x analytically:
    //   μ(x) = c + x*d  where x ~ N(μ_x, σ²_x)
    //   E[y] = c + μ_x * d
    //   Cov(y) = V + σ²_x * d * d'  (rank-1 update)

    for (k in 1:G_mis) {
        int g = index_mis[k];
        int n_g = len[g];
        int start_idx = pos[g];

        // Extract subject's data
        vector[n_g] y_g = segment(y, start_idx, n_g);
        vector[n_g] t_g = segment(t, start_idx, n_g);

        // Compute base mean (without x terms)
        vector[n_g] c_g;
        for (i in 1:n_g) {
            int obs_idx = start_idx + i - 1;
            c_g[i] = alpha_main + beta_t * t_g[i] + dot_product(X[obs_idx], beta);
        }

        // x coefficient vector: d[i] = beta_x + beta_x_t * t[i]
        vector[n_g] d_g = rep_vector(beta_x, n_g) + beta_x_t_interaction * t_g;

        // Marginal mean: E[y] = c + E[x] * d
        real mu_x = x_imputation_mean[g];
        vector[n_g] mu_g = c_g + mu_x * d_g;

        // Base marginal covariance (integrating out RE)
        matrix[n_g, n_g] L_V = compute_marginal_cov_chol(t_g, L_Sigma, sigma_main);

        // Augmented covariance: V_aug = V + σ²_x * d * d' (rank-1 update)
        matrix[n_g, n_g] L_V_aug = chol_rank1_update(L_V, d_g, sigma_imputation);

        // Marginal likelihood (integrating out both RE and x)
        y_g ~ multi_normal_cholesky(mu_g, L_V_aug);
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
        real x_g = x_obs[k];
        int n_g = len[g];
        int start_idx = pos[g];

        vector[n_g] y_g = segment(y, start_idx, n_g);
        vector[n_g] t_g = segment(t, start_idx, n_g);

        vector[n_g] c_g;
        for (i in 1:n_g) {
            int obs_idx = start_idx + i - 1;
            c_g[i] = alpha_main + beta_t * t_g[i] + dot_product(X[obs_idx], beta);
        }
        vector[n_g] d_g = rep_vector(beta_x, n_g) + beta_x_t_interaction * t_g;
        vector[n_g] mu_g = c_g + x_g * d_g;

        matrix[n_g, n_g] L_V = compute_marginal_cov_chol(t_g, L_Sigma, sigma_main);

        real lp_y = multi_normal_cholesky_lpdf(y_g | mu_g, L_V);
        real lp_x = normal_lpdf(x_g | x_imputation_mean[g], sigma_imputation);

        log_lik[g] = lp_y + lp_x;
    }

    // Compute log likelihoods for missing subjects
    for (k in 1:G_mis) {
        int g = index_mis[k];
        int n_g = len[g];
        int start_idx = pos[g];

        vector[n_g] y_g = segment(y, start_idx, n_g);
        vector[n_g] t_g = segment(t, start_idx, n_g);

        vector[n_g] c_g;
        for (i in 1:n_g) {
            int obs_idx = start_idx + i - 1;
            c_g[i] = alpha_main + beta_t * t_g[i] + dot_product(X[obs_idx], beta);
        }
        vector[n_g] d_g = rep_vector(beta_x, n_g) + beta_x_t_interaction * t_g;

        real mu_x = x_imputation_mean[g];
        vector[n_g] mu_g = c_g + mu_x * d_g;

        matrix[n_g, n_g] L_V = compute_marginal_cov_chol(t_g, L_Sigma, sigma_main);
        matrix[n_g, n_g] L_V_aug = chol_rank1_update(L_V, d_g, sigma_imputation);

        // Marginal log likelihood
        log_lik[g] = multi_normal_cholesky_lpdf(y_g | mu_g, L_V_aug);
    }
}
