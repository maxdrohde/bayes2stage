// Mixed Effects Model with Marginalized Random Effects AND Marginalized Beta-Binomial X
//
// This model integrates out BOTH:
//   1. Random effects (analytically via MVN)
//   2. Missing count x values (via sum over Beta-Binomial support)
//
// Result: ZERO latent variables for subjects with missing x
//
// For subjects with missing x:
//   p(y_g | theta) = sum_{k=0}^{m} p(y_g | x=k, theta) * BetaBinomial(k | m, alpha, beta)
//
// where p(y_g | x, theta) is MVN after marginalizing out random effects.
//
// Imputation Model Parameterization (mean-precision):
//   mu_g = inv_logit(gamma_0 + w_g' * gamma)   (mean proportion)
//   phi = precision parameter
//   alpha_g = mu_g * phi
//   beta_g = (1 - mu_g) * phi
//   x_g ~ BetaBinomial(m, alpha_g, beta_g)
//
// Trade-off:
//   - More expensive likelihood (m+1 MVN densities per subject with missing x)
//   - But eliminates ALL latent variables (no funnel, no discrete sampling issues)
//   - Handles overdispersion in count data (unlike Binomial)
//
// Best for: Large G with bounded count x causing mixing problems
//
// Note: For m=1 (binary), this reduces to a Beta-Bernoulli model.
// For large m, consider whether the computational cost (m+1 terms) is acceptable.

functions {
    /**
     * Compute marginal covariance matrix for a subject
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
     * Compute mean vector for subject g given x value
     *
     * @param x Value of count covariate (integer, but passed as real for flexibility)
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

    /**
     * Compute log PMF of Beta-Binomial using recurrence relation
     * This is more efficient than calling beta_binomial_lpmf repeatedly
     *
     * BetaBinomial(k | n, alpha, beta) = C(n,k) * B(k+alpha, n-k+beta) / B(alpha, beta)
     *
     * Recurrence: P(k)/P(k-1) = (n-k+1)(k-1+alpha) / [k(n-k+beta)]
     *
     * In log form:
     *   log P(k) = log P(k-1) + log(n-k+1) + log(k-1+alpha) - log(k) - log(n-k+beta)
     *
     * Base case:
     *   log P(0) = lgamma(n+beta) + lgamma(alpha+beta) - lgamma(n+alpha+beta) - lgamma(beta)
     *
     * @param n Number of trials
     * @param alpha Shape parameter alpha > 0
     * @param beta_param Shape parameter beta > 0
     * @return Vector of log probabilities for k = 0, 1, ..., n
     */
    vector beta_binomial_log_pmf_vec(int n, real alpha, real beta_param) {
        vector[n + 1] log_pmf;

        // Base case: k = 0
        // P(0) = B(alpha, n+beta) / B(alpha, beta)
        //      = Gamma(alpha)Gamma(n+beta)Gamma(alpha+beta) / [Gamma(alpha+n+beta)Gamma(alpha)Gamma(beta)]
        //      = Gamma(n+beta)Gamma(alpha+beta) / [Gamma(n+alpha+beta)Gamma(beta)]
        log_pmf[1] = lgamma(n + beta_param) + lgamma(alpha + beta_param)
                     - lgamma(n + alpha + beta_param) - lgamma(beta_param);

        // Recurrence for k = 1, 2, ..., n
        for (k in 1:n) {
            // P(k)/P(k-1) = (n-k+1)(k-1+alpha) / [k(n-k+beta)]
            log_pmf[k + 1] = log_pmf[k]
                             + log(n - k + 1.0) + log(k - 1.0 + alpha)
                             - log(k * 1.0) - log(n - k + beta_param);
        }

        return log_pmf;
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
    array[G_obs] int<lower=0> x_obs;               // Observed x values (counts)

    // Number of trials for the Beta-Binomial (per-subject array for interface consistency)
    // NOTE: Marginalization requires all subjects have the same n_trials
    array[G] int<lower=1> n_trials;

    array[N] int<lower=1, upper=G> id;      // Subject index per observation

    // Subject start position and length (for ragged array access)
    array[G] int<lower=1> pos;              // Start index for subject g
    array[G] int<lower=1> len;              // Number of obs for subject g
}

transformed data {
    // For marginalization, we need a common n_trials value
    // Validate that all subjects have the same n_trials
    int m_trials = n_trials[1];
    for (g in 2:G) {
        if (n_trials[g] != m_trials) {
            reject("Marginalized model requires all subjects have same n_trials. ",
                   "Subject 1 has ", m_trials, " but subject ", g, " has ", n_trials[g]);
        }
    }

    // Validate that observed x values are within bounds
    for (k in 1:G_obs) {
        if (x_obs[k] > m_trials) {
            reject("x_obs[", k, "] = ", x_obs[k], " exceeds n_trials = ", m_trials);
        }
    }
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

    // Imputation model (Beta-Binomial via mean-precision parameterization)
    real alpha_imputation;                   // Intercept for logit(mu)
    vector[S] gamma;                         // Coefficients for logit(mu)
    real<lower=0> phi_imputation;            // Precision parameter (alpha + beta)

    // NOTE: No x_mis parameters! Count x is marginalized out.
    // NOTE: No z_re parameters! Random effects are marginalized out.
}

transformed parameters {
    // Cholesky factor of RE covariance: L_Sigma = diag(sigma_re) * L_corr
    matrix[2, 2] L_Sigma = diag_pre_multiply(sigma_re, L_re);

    // Imputation model: logit(mu) = linear predictor
    // mu = E[x/m] = alpha_beta / (alpha_beta + beta_beta)
    vector[G] eta_imputation = alpha_imputation + Z * gamma;
    vector[G] mu_imputation = inv_logit(eta_imputation);

    // Beta distribution parameters for each subject
    // Using mean-precision parameterization: alpha = mu * phi, beta = (1-mu) * phi
    vector<lower=0>[G] alpha_beta = mu_imputation * phi_imputation;
    vector<lower=0>[G] beta_beta = (1 - mu_imputation) * phi_imputation;
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
    phi_imputation ~ gamma(2, 0.1);  // Weakly informative prior for precision

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

        // Imputation model likelihood (Beta-Binomial)
        x_g ~ beta_binomial(m_trials, alpha_beta[g], beta_beta[g]);
    }

    // =========================================================================
    // Likelihood for Subjects with MISSING x (Fully Marginalized)
    // =========================================================================
    //
    // For these subjects:
    //   p(y_g | theta) = sum_{k=0}^{m} p(y_g | x=k, theta) * p(x=k | theta)
    //
    // This is a mixture of m+1 MVN distributions.
    // We use log_sum_exp for numerical stability.
    //
    // OPTIMIZATION: Factor out linear algebra from the inner loop.
    // The MVN log-density is: const - 0.5 * ||L^{-1}(y - mu)||^2
    // Since mu_j = base_mu + j*d, we have y - mu_j = r - j*d
    // where r = y - base_mu and d = beta_x + beta_x_t * t
    //
    // Pre-solving r_w = L^{-1}*r and d_w = L^{-1}*d (O(n^2) once),
    // the quadratic form becomes: K1 - 2*j*K2 + j^2*K3 (O(1) per j)
    //
    // Complexity: O(n^2 + m) instead of O(m * n^2)

    for (k in 1:G_mis) {
        int g = index_mis[k];
        int n_g = len[g];
        int start_idx = pos[g];

        // Extract subject's data
        vector[n_g] y_g = segment(y, start_idx, n_g);
        vector[n_g] t_g = segment(t, start_idx, n_g);
        matrix[n_g, P] X_g = X[start_idx:(start_idx + n_g - 1), ];

        // Compute marginal covariance (Cholesky factor)
        matrix[n_g, n_g] L_V = compute_marginal_cov_chol(t_g, L_Sigma, sigma_main);

        // Compute base mean (x = 0) and slope vector
        vector[n_g] base_mu = rep_vector(alpha_main, n_g)
                              + beta_t * t_g
                              + X_g * beta;
        vector[n_g] d = rep_vector(beta_x, n_g) + beta_x_t_interaction * t_g;

        // Compute residual from base mean
        vector[n_g] r = y_g - base_mu;

        // Solve linear systems: r_w = L^{-1} * r, d_w = L^{-1} * d
        vector[n_g] r_w = mdivide_left_tri_low(L_V, r);
        vector[n_g] d_w = mdivide_left_tri_low(L_V, d);

        // Pre-compute sufficient statistics for quadratic form
        real K1 = dot_self(r_w);           // ||r_w||^2
        real K2 = dot_product(r_w, d_w);   // r_w' * d_w
        real K3 = dot_self(d_w);           // ||d_w||^2

        // Log determinant of L_V (for MVN normalization constant)
        real log_det_L = sum(log(diagonal(L_V)));

        // MVN normalization constant: -n/2 * log(2*pi) - log|L|
        real mvn_const = -0.5 * n_g * log(2 * pi()) - log_det_L;

        // Compute all Beta-Binomial log PMFs using efficient recurrence
        vector[m_trials + 1] log_pmf_bb = beta_binomial_log_pmf_vec(
            m_trials, alpha_beta[g], beta_beta[g]);

        // Compute log probability for each possible x value
        vector[m_trials + 1] log_probs;

        for (j in 0:m_trials) {
            // Quadratic form: ||L^{-1}(y - mu_j)||^2 = K1 - 2*j*K2 + j^2*K3
            real Q = K1 - 2.0 * j * K2 + square(j) * K3;

            // Log density of outcome given x = j
            real lp_y_xj = mvn_const - 0.5 * Q;

            // Log probability of x = j from Beta-Binomial (pre-computed)
            log_probs[j + 1] = lp_y_xj + log_pmf_bb[j + 1];
        }

        // Marginal likelihood via log-sum-exp
        target += log_sum_exp(log_probs);
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
        real lp_x = beta_binomial_lpmf(x_g | m_trials, alpha_beta[g], beta_beta[g]);

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

        // Compute base mean and slope vector
        vector[n_g] base_mu = rep_vector(alpha_main, n_g)
                              + beta_t * t_g
                              + X_g * beta;
        vector[n_g] d = rep_vector(beta_x, n_g) + beta_x_t_interaction * t_g;
        vector[n_g] r = y_g - base_mu;

        // Solve linear systems once
        vector[n_g] r_w = mdivide_left_tri_low(L_V, r);
        vector[n_g] d_w = mdivide_left_tri_low(L_V, d);

        // Pre-compute sufficient statistics
        real K1 = dot_self(r_w);
        real K2 = dot_product(r_w, d_w);
        real K3 = dot_self(d_w);
        real log_det_L = sum(log(diagonal(L_V)));
        real mvn_const = -0.5 * n_g * log(2 * pi()) - log_det_L;

        // Compute all Beta-Binomial log PMFs using efficient recurrence
        vector[m_trials + 1] log_pmf_bb = beta_binomial_log_pmf_vec(
            m_trials, alpha_beta[g], beta_beta[g]);

        vector[m_trials + 1] log_probs;

        for (j in 0:m_trials) {
            real Q = K1 - 2.0 * j * K2 + square(j) * K3;
            real lp_y_xj = mvn_const - 0.5 * Q;
            log_probs[j + 1] = lp_y_xj + log_pmf_bb[j + 1];
        }

        // Marginal log likelihood
        log_lik[g] = log_sum_exp(log_probs);
    }
}
