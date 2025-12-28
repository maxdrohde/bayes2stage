// Double-Woodbury optimized Stan model
// Diagonal + low-rank covariance throughout:
//   V = sigma^2 I + Z * Sigma_re * Z'                (rank 2)
// Missing x adds: + sigma_x^2 * d d'                 (rank 1)
// => observed-x subjects: diagonal + rank-2
// => missing-x subjects:  diagonal + rank-3 (per mixture component)
//
// Key optimization: never form n×n V and never Cholesky it.
// We only Cholesky tiny k×k matrices (k=2 or 3).
//
// Complexity:
//   Current model: O(n³) per subject for Cholesky
//   This model:    O(n) per subject for sufficient stats + O(k³) for k×k Cholesky
//
// Best for: Large n_g (many observations per subject)

functions {
    /**
     * MVN log density for covariance: Sigma = sigma_sq * I + Q Q'
     * using only sufficient statistics:
     *   r2  = r'r
     *   QQ  = Q'Q  (k×k)
     *   Qtr = Q'r  (k)
     *
     * Uses matrix determinant lemma and Woodbury identity:
     *   log|Sigma| = n*log(sigma_sq) + log|I + Q'Q/sigma_sq|
     *   r'Sigma^{-1}r = r'r/sigma_sq - (Q'r)'M^{-1}(Q'r)/sigma_sq^2
     *   where M = I + Q'Q/sigma_sq
     *
     * Complexity: O(k^3) for Cholesky of k×k, with k=2 or 3.
     *
     * @param r2 Scalar r'r (sum of squared residuals)
     * @param Qtr Vector Q'r (k-vector)
     * @param QQ Matrix Q'Q (k×k)
     * @param sigma_sq Residual variance
     * @param n Number of observations
     * @return Log probability density
     */
    real mvn_diag_plus_lowrank_lpdf(real r2,
                                    vector Qtr,
                                    matrix QQ,
                                    real sigma_sq,
                                    int n) {
        int k = rows(Qtr);

        // M = I + Q'Q / sigma_sq
        matrix[k, k] M = QQ / sigma_sq;
        for (i in 1:k) {
            M[i, i] += 1.0;
        }

        matrix[k, k] L_M = cholesky_decompose(M);

        // log|sigma_sq * I + Q Q'| = n * log(sigma_sq) + log|M|
        //                          = n * log(sigma_sq) + 2 * sum(log(diag(L_M)))
        real log_det = n * log(sigma_sq) + 2.0 * sum(log(diagonal(L_M)));

        // Quadratic form: r' Sigma^{-1} r
        //   = r'r / sigma_sq - (Q'r)' M^{-1} (Q'r) / sigma_sq^2
        vector[k] sol = mdivide_left_tri_low(L_M, Qtr);
        real quad = (r2 - dot_self(sol) / sigma_sq) / sigma_sq;

        return -0.5 * (n * log(2.0 * pi()) + log_det + quad);
    }
}

data {
    // Dimensions
    int<lower=1> N;                         // Total observations
    int<lower=1> G;                         // Total subjects
    int<lower=0, upper=G> G_obs;            // Subjects with observed x
    int<lower=0, upper=G> G_mis;            // Subjects with missing x

    int<lower=0> P;                         // Main model predictors (excluding x)
    vector[N] t;                            // Time variable
    matrix[N, P] X;                         // Main model covariates
    vector[N] y;                            // Outcome variable

    int<lower=0> S;                         // Imputation model predictors
    matrix[G, S] Z;                         // Imputation covariates (subject-level)

    // Subject indexing
    array[G_obs] int<lower=1, upper=G> index_obs;  // Subjects with observed x
    array[G_mis] int<lower=1, upper=G> index_mis;  // Subjects with missing x
    vector[G_obs] x_obs;                           // Observed x values

    // Ragged array access
    array[G] int<lower=1> pos;              // Start index for subject g
    array[G] int<lower=1> len;              // Number of obs for subject g

    // Mixture model settings
    int<lower=2> K;                         // Number of mixture components
}

parameters {
    // Random effects covariance (intercept + slope)
    vector<lower=0>[2] sigma_re;            // [intercept SD, slope SD]
    cholesky_factor_corr[2] L_re;           // Cholesky factor of RE correlation

    // Main model fixed effects
    real alpha_main;
    real beta_t;
    real beta_x;
    real beta_x_t_interaction;
    vector[P] beta;
    real<lower=1e-8> sigma_main;

    // Mixture imputation model for x
    simplex[K] theta;
    real alpha_imputation;
    vector[S] gamma;
    ordered[K] mu_component_raw;            // Ordered offsets (will be centered)
    vector<lower=0>[K] sigma_component;
}

transformed parameters {
    // Lower Cholesky of RE covariance: L_Sigma = diag(sigma_re) * L_re
    matrix[2, 2] L_Sigma = diag_pre_multiply(sigma_re, L_re);

    // Precomputed quantities used in model block
    real sigma_sq = square(sigma_main);
    vector[N] xb = X * beta;

    // Center component offsets for identifiability
    vector[K] mu_component = mu_component_raw - mean(mu_component_raw);

    // Subject baseline mean for x
    vector[G] x_base_mean = alpha_imputation + Z * gamma;

    // L_Sigma elements for sufficient stats computation
    // L_Sigma = [[a, 0], [b, c]]
    real re_a = L_Sigma[1, 1];
    real re_b = L_Sigma[2, 1];
    real re_c = L_Sigma[2, 2];
}

model {
    // =========================================================================
    // Priors
    // =========================================================================
    sigma_re ~ exponential(0.1);
    L_re ~ lkj_corr_cholesky(2);

    sigma_main ~ exponential(0.1);
    alpha_main ~ normal(0, 100);
    beta_t ~ normal(0, 100);
    beta_x ~ normal(0, 100);
    beta_x_t_interaction ~ normal(0, 100);
    beta ~ normal(0, 100);

    theta ~ dirichlet(rep_vector(2.0, K));
    alpha_imputation ~ normal(0, 100);
    gamma ~ normal(0, 100);
    mu_component_raw ~ normal(0, 10);
    sigma_component ~ exponential(0.1);

    // =========================================================================
    // OBSERVED x subjects: diagonal + rank-2 (random effects only)
    // =========================================================================
    for (k_subj in 1:G_obs) {
        int g = index_obs[k_subj];
        int n_g = len[g];
        int start_idx = pos[g];
        real x_g = x_obs[k_subj];

        // Accumulate sufficient statistics in one pass: O(n_g)
        real sumt = 0;
        real sumt2 = 0;
        real sum_r = 0;
        real sum_tr = 0;
        real sum_r2 = 0;

        for (i in 1:n_g) {
            int obs = start_idx + i - 1;
            real ti = t[obs];
            real mu_base_i = alpha_main + beta_t * ti + xb[obs];
            real ri = y[obs] - mu_base_i;

            sumt += ti;
            sumt2 += ti * ti;
            sum_r += ri;
            sum_tr += ti * ri;
            sum_r2 += ri * ri;
        }

        // d(t) = beta_x + beta_x_t_interaction * t
        // d'd = sum_i d(t_i)^2
        real dd = n_g * square(beta_x)
                  + 2.0 * beta_x * beta_x_t_interaction * sumt
                  + square(beta_x_t_interaction) * sumt2;

        // d' r_base
        real d_rbase = beta_x * sum_r + beta_x_t_interaction * sum_tr;

        // r = r_base - x * d  =>  r'r
        real r2 = sum_r2 - 2.0 * x_g * d_rbase + square(x_g) * dd;

        // Q_re = Z_re * L_Sigma where Z_re = [1, t]
        // Columns: q1 = a + b*t, q2 = c*t
        // Q_re' Q_re (2×2)
        matrix[2, 2] QQ_re;
        QQ_re[1, 1] = n_g * square(re_a) + 2.0 * re_a * re_b * sumt + square(re_b) * sumt2;
        QQ_re[1, 2] = re_c * re_a * sumt + re_c * re_b * sumt2;
        QQ_re[2, 1] = QQ_re[1, 2];
        QQ_re[2, 2] = square(re_c) * sumt2;

        // Q_re' d (used for Q_re' r = Q_re' r_base - x * Q_re' d)
        vector[2] Qd;
        Qd[1] = n_g * re_a * beta_x
                + (re_a * beta_x_t_interaction + re_b * beta_x) * sumt
                + re_b * beta_x_t_interaction * sumt2;
        Qd[2] = re_c * beta_x * sumt + re_c * beta_x_t_interaction * sumt2;

        // Q_re' r_base
        vector[2] Qtr_base;
        Qtr_base[1] = re_a * sum_r + re_b * sum_tr;
        Qtr_base[2] = re_c * sum_tr;

        // Q_re' r = Q_re' (r_base - x * d) = Qtr_base - x * Qd
        vector[2] Qtr = Qtr_base - x_g * Qd;

        // Likelihood for y | x observed: O(1) beyond the sums
        target += mvn_diag_plus_lowrank_lpdf(r2 | Qtr, QQ_re, sigma_sq, n_g);

        // Likelihood for observed x under mixture
        vector[K] lps_x;
        for (k in 1:K) {
            real mu_gk = x_base_mean[g] + mu_component[k];
            lps_x[k] = log(theta[k]) + normal_lpdf(x_g | mu_gk, sigma_component[k]);
        }
        target += log_sum_exp(lps_x);
    }

    // =========================================================================
    // MISSING x subjects: diagonal + rank-3 (random effects + missing-x)
    // =========================================================================
    for (k_subj in 1:G_mis) {
        int g = index_mis[k_subj];
        int n_g = len[g];
        int start_idx = pos[g];

        // Accumulate sufficient statistics: O(n_g)
        real sumt = 0;
        real sumt2 = 0;
        real sum_r = 0;
        real sum_tr = 0;
        real sum_r2 = 0;

        for (i in 1:n_g) {
            int obs = start_idx + i - 1;
            real ti = t[obs];
            real mu_base_i = alpha_main + beta_t * ti + xb[obs];
            real ri = y[obs] - mu_base_i;

            sumt += ti;
            sumt2 += ti * ti;
            sum_r += ri;
            sum_tr += ti * ri;
            sum_r2 += ri * ri;
        }

        // Shared quantities for this subject
        real dd = n_g * square(beta_x)
                  + 2.0 * beta_x * beta_x_t_interaction * sumt
                  + square(beta_x_t_interaction) * sumt2;

        real d_rbase = beta_x * sum_r + beta_x_t_interaction * sum_tr;

        // Random effects low-rank pieces (constant across mixture components)
        matrix[2, 2] QQ_re;
        QQ_re[1, 1] = n_g * square(re_a) + 2.0 * re_a * re_b * sumt + square(re_b) * sumt2;
        QQ_re[1, 2] = re_c * re_a * sumt + re_c * re_b * sumt2;
        QQ_re[2, 1] = QQ_re[1, 2];
        QQ_re[2, 2] = square(re_c) * sumt2;

        vector[2] Qd;
        Qd[1] = n_g * re_a * beta_x
                + (re_a * beta_x_t_interaction + re_b * beta_x) * sumt
                + re_b * beta_x_t_interaction * sumt2;
        Qd[2] = re_c * beta_x * sumt + re_c * beta_x_t_interaction * sumt2;

        vector[2] Qtr_base;
        Qtr_base[1] = re_a * sum_r + re_b * sum_tr;
        Qtr_base[2] = re_c * sum_tr;

        // Marginalize over mixture components
        vector[K] lps;
        for (k in 1:K) {
            real mu_k = x_base_mean[g] + mu_component[k];
            real sx = sigma_component[k];

            // r = r_base - mu_k * d  =>  r'r
            real r2 = sum_r2 - 2.0 * mu_k * d_rbase + square(mu_k) * dd;

            // Q_re' r = Qtr_base - mu_k * Qd
            vector[2] Qtr_re = Qtr_base - mu_k * Qd;

            // d' r = d' r_base - mu_k * d'd
            real d_r = d_rbase - mu_k * dd;

            // Build rank-3 sufficient stats:
            // Q_full = [Q_re, sx * d]
            // Q_full' Q_full (3×3):
            matrix[3, 3] QQ_full;
            QQ_full[1:2, 1:2] = QQ_re;
            QQ_full[1, 3] = sx * Qd[1];
            QQ_full[2, 3] = sx * Qd[2];
            QQ_full[3, 1] = sx * Qd[1];
            QQ_full[3, 2] = sx * Qd[2];
            QQ_full[3, 3] = square(sx) * dd;

            // Q_full' r (3-vector):
            vector[3] Qtr_full;
            Qtr_full[1:2] = Qtr_re;
            Qtr_full[3] = sx * d_r;

            lps[k] = log(theta[k]) +
                     mvn_diag_plus_lowrank_lpdf(r2 | Qtr_full, QQ_full, sigma_sq, n_g);
        }

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

    // Observed x subjects
    for (k_subj in 1:G_obs) {
        int g = index_obs[k_subj];
        int n_g = len[g];
        int start_idx = pos[g];
        real x_g = x_obs[k_subj];

        real sumt = 0;
        real sumt2 = 0;
        real sum_r = 0;
        real sum_tr = 0;
        real sum_r2 = 0;

        for (i in 1:n_g) {
            int obs = start_idx + i - 1;
            real ti = t[obs];
            real mu_base_i = alpha_main + beta_t * ti + xb[obs];
            real ri = y[obs] - mu_base_i;

            sumt += ti;
            sumt2 += ti * ti;
            sum_r += ri;
            sum_tr += ti * ri;
            sum_r2 += ri * ri;
        }

        real dd = n_g * square(beta_x)
                  + 2.0 * beta_x * beta_x_t_interaction * sumt
                  + square(beta_x_t_interaction) * sumt2;

        real d_rbase = beta_x * sum_r + beta_x_t_interaction * sum_tr;
        real r2 = sum_r2 - 2.0 * x_g * d_rbase + square(x_g) * dd;

        matrix[2, 2] QQ_re;
        QQ_re[1, 1] = n_g * square(re_a) + 2.0 * re_a * re_b * sumt + square(re_b) * sumt2;
        QQ_re[1, 2] = re_c * re_a * sumt + re_c * re_b * sumt2;
        QQ_re[2, 1] = QQ_re[1, 2];
        QQ_re[2, 2] = square(re_c) * sumt2;

        vector[2] Qd;
        Qd[1] = n_g * re_a * beta_x
                + (re_a * beta_x_t_interaction + re_b * beta_x) * sumt
                + re_b * beta_x_t_interaction * sumt2;
        Qd[2] = re_c * beta_x * sumt + re_c * beta_x_t_interaction * sumt2;

        vector[2] Qtr_base;
        Qtr_base[1] = re_a * sum_r + re_b * sum_tr;
        Qtr_base[2] = re_c * sum_tr;

        vector[2] Qtr = Qtr_base - x_g * Qd;

        real lp_y = mvn_diag_plus_lowrank_lpdf(r2 | Qtr, QQ_re, sigma_sq, n_g);

        vector[K] lps_x;
        for (k in 1:K) {
            real mu_gk = x_base_mean[g] + mu_component[k];
            lps_x[k] = log(theta[k]) + normal_lpdf(x_g | mu_gk, sigma_component[k]);
        }
        real lp_x = log_sum_exp(lps_x);

        log_lik[g] = lp_y + lp_x;
    }

    // Missing x subjects
    for (k_subj in 1:G_mis) {
        int g = index_mis[k_subj];
        int n_g = len[g];
        int start_idx = pos[g];

        real sumt = 0;
        real sumt2 = 0;
        real sum_r = 0;
        real sum_tr = 0;
        real sum_r2 = 0;

        for (i in 1:n_g) {
            int obs = start_idx + i - 1;
            real ti = t[obs];
            real mu_base_i = alpha_main + beta_t * ti + xb[obs];
            real ri = y[obs] - mu_base_i;

            sumt += ti;
            sumt2 += ti * ti;
            sum_r += ri;
            sum_tr += ti * ri;
            sum_r2 += ri * ri;
        }

        real dd = n_g * square(beta_x)
                  + 2.0 * beta_x * beta_x_t_interaction * sumt
                  + square(beta_x_t_interaction) * sumt2;

        real d_rbase = beta_x * sum_r + beta_x_t_interaction * sum_tr;

        matrix[2, 2] QQ_re;
        QQ_re[1, 1] = n_g * square(re_a) + 2.0 * re_a * re_b * sumt + square(re_b) * sumt2;
        QQ_re[1, 2] = re_c * re_a * sumt + re_c * re_b * sumt2;
        QQ_re[2, 1] = QQ_re[1, 2];
        QQ_re[2, 2] = square(re_c) * sumt2;

        vector[2] Qd;
        Qd[1] = n_g * re_a * beta_x
                + (re_a * beta_x_t_interaction + re_b * beta_x) * sumt
                + re_b * beta_x_t_interaction * sumt2;
        Qd[2] = re_c * beta_x * sumt + re_c * beta_x_t_interaction * sumt2;

        vector[2] Qtr_base;
        Qtr_base[1] = re_a * sum_r + re_b * sum_tr;
        Qtr_base[2] = re_c * sum_tr;

        vector[K] lps;
        for (k in 1:K) {
            real mu_k = x_base_mean[g] + mu_component[k];
            real sx = sigma_component[k];

            real r2 = sum_r2 - 2.0 * mu_k * d_rbase + square(mu_k) * dd;

            vector[2] Qtr_re = Qtr_base - mu_k * Qd;
            real d_r = d_rbase - mu_k * dd;

            matrix[3, 3] QQ_full;
            QQ_full[1:2, 1:2] = QQ_re;
            QQ_full[1, 3] = sx * Qd[1];
            QQ_full[2, 3] = sx * Qd[2];
            QQ_full[3, 1] = sx * Qd[1];
            QQ_full[3, 2] = sx * Qd[2];
            QQ_full[3, 3] = square(sx) * dd;

            vector[3] Qtr_full;
            Qtr_full[1:2] = Qtr_re;
            Qtr_full[3] = sx * d_r;

            lps[k] = log(theta[k]) +
                     mvn_diag_plus_lowrank_lpdf(r2 | Qtr_full, QQ_full, sigma_sq, n_g);
        }

        log_lik[g] = log_sum_exp(lps);
    }
}

