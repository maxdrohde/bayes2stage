// Mixed Effects Model with Fully Marginalized Random Effects AND Missing X
// with MIXTURE OF NORMALS for the imputation model
//
// This model integrates out BOTH:
//   1. Random effects (analytically via MVN)
//   2. Missing continuous x values (analytically via Gaussian mixture)
//
// Result: ZERO latent variables
//
// For subjects with missing x and mixture model:
//   x ~ Σ_k π_k * N(μ_k, σ²_k)
//   p(y | x missing) = Σ_k π_k * p(y | x ~ N(μ_k, σ²_k))
//
// Each component contributes a marginalized MVN:
//   y | component k ~ MVN(c + μ_k*d, V + σ²_k*d*d')
//
// Trade-off:
//   - More flexible imputation model for multimodal x distributions
//   - Computation scales with number of mixture components
//   - Still eliminates ALL latent variables
//
// Best for: Large G with continuous x that has multimodal distribution

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
     * MVN log-density with rank-1 covariance update using Woodbury identity
     *
     * Computes log p(y | mu, V + σ² d d') efficiently using:
     *   - Determinant lemma: log|V + σ²dd'| = log|V| + log(1 + σ² d'V⁻¹d)
     *   - Woodbury identity: (V + σ²dd')⁻¹ = V⁻¹ - σ²/(1 + σ²d'V⁻¹d) (V⁻¹d)(V⁻¹d)'
     *
     * Truly O(n) per component when all precomputation is done:
     *   - z = L\d, w_base = L\(y-c), z_norm_sq, log_det_V computed once
     *   - For component k: w_k = w_base - μ_k*z is O(n), then dot products O(n)
     *
     * @param w Precomputed L \ r where r = y - mu (length n)
     * @param z Precomputed L \ d (length n)
     * @param z_norm_sq Precomputed ||z||² = d' V⁻¹ d
     * @param log_det_V Precomputed log|V| = 2 * sum(log(diag(L)))
     * @param sigma_x Standard deviation of x for rank-1 update
     * @param n Dimension
     * @return Log probability density
     */
    real mvn_rank1_woodbury_lpdf(vector w, vector z, real z_norm_sq,
                                  real log_det_V, real sigma_x, int n) {
        real sigma_sq = square(sigma_x);

        // Woodbury denominator: κ = 1 + σ² ||z||²
        real kappa = 1 + sigma_sq * z_norm_sq;

        // Log determinant: log|V_aug| = log|V| + log(κ)
        real log_det_V_aug = log_det_V + log(kappa);

        // Quadratic form using Woodbury:
        // r' V_aug⁻¹ r = r' V⁻¹ r - σ²/κ (r' V⁻¹ d)²
        //              = ||w||² - σ²/κ (w' z)²
        real w_norm_sq = dot_self(w);
        real wz = dot_product(w, z);
        real quad_form = w_norm_sq - sigma_sq / kappa * square(wz);

        // MVN log-density: -n/2 log(2π) - 1/2 log|V_aug| - 1/2 quad
        return -0.5 * n * log(2 * pi()) - 0.5 * log_det_V_aug - 0.5 * quad_form;
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

    // Mixture model settings
    int<lower=2> K;                         // Number of mixture components
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

    // Mixture imputation model
    simplex[K] theta;                       // Mixture proportions
    real alpha_imputation;                  // Intercept (overall mean of mixture)
    vector[S] gamma;                        // Covariate effects (shared across components)
    ordered[K] mu_component_raw;            // Ordered component offsets (unconstrained)
    vector<lower=0>[K] sigma_component;     // Component-specific standard deviations

    // NOTE: No x_mis parameters! Continuous x is marginalized out.
    // NOTE: No z_re parameters! Random effects are marginalized out.
}

transformed parameters {
    // Cholesky factor of RE covariance: L_Sigma = diag(sigma_re) * L_corr
    matrix[2, 2] L_Sigma = diag_pre_multiply(sigma_re, L_re);

    // Center mu_component at zero while preserving ordering
    // This ensures identifiability: alpha_imputation captures overall mean,
    // mu_component captures deviations with sum(mu_component) = 0
    // Ordering is preserved since we subtract a constant from all elements
    vector[K] mu_component = mu_component_raw - mean(mu_component_raw);

    // Imputation model: subject-specific base mean (before component offset)
    // For component k: mean = x_base_mean[g] + mu_component[k]
    vector[G] x_base_mean = alpha_imputation + Z * gamma;
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

    // Mixture imputation model priors
    theta ~ dirichlet(rep_vector(2.0, K));  // Slightly informative prior toward equal weights
    alpha_imputation ~ normal(0, 100);      // Overall mean of mixture
    gamma ~ normal(0, 100);
    mu_component_raw ~ normal(0, 10);       // Component offsets (ordering enforced by constraint)
    sigma_component ~ exponential(0.1);

    // =========================================================================
    // Likelihood for Subjects with OBSERVED x
    // =========================================================================

    for (k_subj in 1:G_obs) {
        int g = index_obs[k_subj];
        real x_g = x_obs[k_subj];
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

        // Mixture imputation model likelihood for observed x
        // p(x_g) = Σ_k θ_k * N(x_g | μ_gk, σ_k)
        vector[K] lps;
        for (k_comp in 1:K) {
            real mu_gk = x_base_mean[g] + mu_component[k_comp];
            lps[k_comp] = log(theta[k_comp]) + normal_lpdf(x_g | mu_gk, sigma_component[k_comp]);
        }
        target += log_sum_exp(lps);
    }

    // =========================================================================
    // Likelihood for Subjects with MISSING x (Fully Marginalized)
    // =========================================================================
    //
    // For these subjects with mixture model on x:
    //   p(y | x missing) = Σ_k θ_k * p(y | x ~ N(μ_gk, σ²_k))
    //
    // Each component contributes:
    //   E[y | component k] = c + μ_gk * d
    //   Cov(y | component k) = V + σ²_k * d * d'
    //
    // Uses Woodbury identity for O(n) per component instead of O(n²) Cholesky update.

    for (k_subj in 1:G_mis) {
        int g = index_mis[k_subj];
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

        // Base marginal covariance (integrating out RE)
        matrix[n_g, n_g] L_V = compute_marginal_cov_chol(t_g, L_Sigma, sigma_main);

        // Precompute quantities for Woodbury identity (O(n²) once, then O(n) per component)
        vector[n_g] z = mdivide_left_tri_low(L_V, d_g);       // z = L \ d
        real z_norm_sq = dot_self(z);                         // ||z||² = d' V⁻¹ d
        real log_det_V = 2 * sum(log(diagonal(L_V)));         // log|V|
        vector[n_g] r_base = y_g - c_g;                       // Base residual
        vector[n_g] w_base = mdivide_left_tri_low(L_V, r_base); // w_base = L \ (y - c)

        // Marginalize over mixture components using Woodbury (truly O(n) per component)
        // Key: L \ (r_base - μ*d) = w_base - μ*z by linearity
        vector[K] lps;
        for (k_comp in 1:K) {
            real mu_gk = x_base_mean[g] + mu_component[k_comp];

            // w_k = L \ r_k = L \ (r_base - μ_gk * d) = w_base - μ_gk * z  (O(n))
            vector[n_g] w_k = w_base - mu_gk * z;

            // Log-probability using Woodbury identity (O(n) - just dot products)
            lps[k_comp] = log(theta[k_comp]) +
                          mvn_rank1_woodbury_lpdf(w_k | z, z_norm_sq, log_det_V,
                                                   sigma_component[k_comp], n_g);
        }

        // Marginal likelihood (log-sum-exp over components)
        target += log_sum_exp(lps);
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
    for (k_subj in 1:G_obs) {
        int g = index_obs[k_subj];
        real x_g = x_obs[k_subj];
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

        // Mixture log-likelihood for x
        vector[K] lps_x;
        for (k_comp in 1:K) {
            real mu_gk = x_base_mean[g] + mu_component[k_comp];
            lps_x[k_comp] = log(theta[k_comp]) + normal_lpdf(x_g | mu_gk, sigma_component[k_comp]);
        }
        real lp_x = log_sum_exp(lps_x);

        log_lik[g] = lp_y + lp_x;
    }

    // Compute log likelihoods for missing subjects (using Woodbury identity)
    for (k_subj in 1:G_mis) {
        int g = index_mis[k_subj];
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

        matrix[n_g, n_g] L_V = compute_marginal_cov_chol(t_g, L_Sigma, sigma_main);

        // Precompute for Woodbury
        vector[n_g] z = mdivide_left_tri_low(L_V, d_g);
        real z_norm_sq = dot_self(z);
        real log_det_V = 2 * sum(log(diagonal(L_V)));
        vector[n_g] r_base = y_g - c_g;
        vector[n_g] w_base = mdivide_left_tri_low(L_V, r_base);

        vector[K] lps;
        for (k_comp in 1:K) {
            real mu_gk = x_base_mean[g] + mu_component[k_comp];
            vector[n_g] w_k = w_base - mu_gk * z;
            lps[k_comp] = log(theta[k_comp]) +
                          mvn_rank1_woodbury_lpdf(w_k | z, z_norm_sq, log_det_V,
                                                   sigma_component[k_comp], n_g);
        }

        log_lik[g] = log_sum_exp(lps);
    }
}

