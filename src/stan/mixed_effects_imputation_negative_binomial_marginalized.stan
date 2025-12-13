// Mixed Effects Model with Marginalized Random Effects AND Marginalized Negative Binomial X
//
// This model integrates out BOTH:
//   1. Random effects (analytically via MVN)
//   2. Missing count x values (via truncated sum over Negative Binomial support)
//
// Result: ZERO latent variables for subjects with missing x
//
// CRITICAL DIFFERENCE from Beta-Binomial:
//   Negative Binomial has UNBOUNDED support [0, infinity).
//   We must TRUNCATE the sum at x_max, introducing approximation error.
//
// For subjects with missing x:
//   p(y_g | theta) ≈ sum_{k=0}^{x_max} p(y_g | x=k, theta) * NegBin(k | mu, phi)
//
// where p(y_g | x, theta) is MVN after marginalizing out random effects.
//
// Imputation Model Parameterization (log-link, mean-overdispersion):
//   log(mu_g) = gamma_0 + w_g' * gamma
//   x_g ~ NegBinomial2(mu_g, phi)
//   E[x] = mu, Var[x] = mu + mu^2/phi
//
// Truncation considerations:
//   - x_max should be large enough that P(x > x_max) is negligible
//   - Rule of thumb: x_max = max(observed_x) + buffer, or based on domain knowledge
//   - The model outputs truncation diagnostics in generated quantities
//
// Trade-off:
//   - More expensive likelihood (x_max+1 MVN densities per subject with missing x)
//   - Truncation introduces small approximation error
//   - But eliminates ALL latent variables (no funnel, no discrete sampling issues)
//
// Best for: Large G with overdispersed count x causing mixing problems

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
     * @param x Value of count covariate
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
     * Compute log PMF of Negative Binomial using recurrence relation
     * This is more numerically stable than calling neg_binomial_2_lpmf repeatedly
     *
     * NegBin2(x | mu, phi) parameterization:
     *   E[X] = mu
     *   Var[X] = mu + mu^2/phi
     *
     * Recurrence: P(k) = P(k-1) * (k-1+phi)/k * mu/(mu+phi)
     *
     * @param x_max Maximum x value to compute
     * @param mu Mean parameter
     * @param phi Overdispersion parameter
     * @return Vector of log probabilities for x = 0, 1, ..., x_max
     */
    vector negbinomial_log_pmf_vec(int x_max, real mu, real phi) {
        vector[x_max + 1] log_pmf;

        // Base case: x = 0
        // P(0) = (phi / (mu + phi))^phi
        log_pmf[1] = phi * log(phi / (mu + phi));

        // Recurrence for x = 1, 2, ..., x_max
        real log_ratio = log(mu / (mu + phi));
        for (k in 1:x_max) {
            log_pmf[k + 1] = log_pmf[k] + log(k - 1.0 + phi) - log(k) + log_ratio;
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

    // Truncation for Negative Binomial marginalization
    // IMPORTANT: Must be >= max(x_obs) and large enough for good approximation
    int<lower=0> x_max;

    array[N] int<lower=1, upper=G> id;      // Subject index per observation

    // Subject start position and length (for ragged array access)
    array[G] int<lower=1> pos;              // Start index for subject g
    array[G] int<lower=1> len;              // Number of obs for subject g
}

transformed data {
    // Validate that x_max is at least as large as observed values
    int max_x_obs = 0;
    for (k in 1:G_obs) {
        if (x_obs[k] > max_x_obs) {
            max_x_obs = x_obs[k];
        }
    }
    if (x_max < max_x_obs) {
        reject("x_max (", x_max, ") must be >= max(x_obs) (", max_x_obs, ")");
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

    // Imputation model (Negative Binomial via log-link)
    real alpha_imputation;                   // Intercept for log(mu)
    vector[S] gamma;                         // Coefficients for log(mu)
    real<lower=0> phi_imputation;            // Overdispersion (larger = less overdispersion)

    // NOTE: No x_mis parameters! Count x is marginalized out.
    // NOTE: No z_re parameters! Random effects are marginalized out.
}

transformed parameters {
    // Cholesky factor of RE covariance: L_Sigma = diag(sigma_re) * L_corr
    matrix[2, 2] L_Sigma = diag_pre_multiply(sigma_re, L_re);

    // Imputation model: log(mu) = linear predictor
    vector[G] eta_imputation = alpha_imputation + Z * gamma;
    vector<lower=0>[G] mu_imputation = exp(eta_imputation);
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
    phi_imputation ~ gamma(2, 0.1);  // Weakly informative prior for overdispersion

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

        // Imputation model likelihood (Negative Binomial)
        x_g ~ neg_binomial_2(mu_imputation[g], phi_imputation);
    }

    // =========================================================================
    // Likelihood for Subjects with MISSING x (Truncated Marginalization)
    // =========================================================================
    //
    // For these subjects:
    //   p(y_g | theta) ≈ sum_{k=0}^{x_max} p(y_g | x=k, theta) * p(x=k | theta)
    //
    // This is a mixture of x_max+1 MVN distributions.
    // We use log_sum_exp for numerical stability.
    //
    // NOTE: This is an APPROXIMATION due to truncation.
    // The error is P(x > x_max), which should be negligible if x_max is large enough.
    //
    // OPTIMIZATION: Factor out linear algebra from the inner loop.
    // The MVN log-density is: const - 0.5 * ||L^{-1}(y - mu)||^2
    // Since mu_j = base_mu + j*d, we have y - mu_j = r - j*d
    // where r = y - base_mu and d = beta_x + beta_x_t * t
    //
    // Pre-solving r_w = L^{-1}*r and d_w = L^{-1}*d (O(n^2) once),
    // the quadratic form becomes: K1 - 2*j*K2 + j^2*K3 (O(1) per j)
    //
    // Complexity: O(n^2 + x_max) instead of O(x_max * n^2)

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

        // Precompute NegBin log PMF using recurrence (more efficient)
        vector[x_max + 1] log_pmf_x = negbinomial_log_pmf_vec(x_max,
                                                              mu_imputation[g],
                                                              phi_imputation);

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

        // Compute log probability for each possible x value
        vector[x_max + 1] log_probs;

        for (j in 0:x_max) {
            // Quadratic form: ||L^{-1}(y - mu_j)||^2 = K1 - 2*j*K2 + j^2*K3
            real Q = K1 - 2.0 * j * K2 + square(j) * K3;

            // Log density of outcome given x = j
            real lp_y_xj = mvn_const - 0.5 * Q;

            // Log probability of x = j from Negative Binomial imputation model
            real lp_xj = log_pmf_x[j + 1];

            log_probs[j + 1] = lp_y_xj + lp_xj;
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

    // Truncation diagnostics: P(x > x_max) for each missing subject
    // If these are not small, x_max should be increased
    vector[G_mis] prob_x_exceeds_max;

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
        real lp_x = neg_binomial_2_lpmf(x_g | mu_imputation[g], phi_imputation);

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

        // Precompute NegBin log PMF
        vector[x_max + 1] log_pmf_x = negbinomial_log_pmf_vec(x_max,
                                                              mu_imputation[g],
                                                              phi_imputation);

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

        vector[x_max + 1] log_probs;

        for (j in 0:x_max) {
            real Q = K1 - 2.0 * j * K2 + square(j) * K3;
            real lp_y_xj = mvn_const - 0.5 * Q;
            log_probs[j + 1] = lp_y_xj + log_pmf_x[j + 1];
        }

        // Marginal log likelihood
        log_lik[g] = log_sum_exp(log_probs);

        // Truncation diagnostic: P(x > x_max | imputation model only)
        // This is an upper bound on the posterior truncation error
        // Compute as 1 - sum_{k=0}^{x_max} P(X = k)
        real cumulative_prob = exp(log_sum_exp(log_pmf_x));
        prob_x_exceeds_max[k] = 1.0 - cumulative_prob;
    }
}
