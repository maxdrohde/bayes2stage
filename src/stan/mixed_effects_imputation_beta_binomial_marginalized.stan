// Mixed Effects Model with Marginalized Random Effects AND Marginalized Beta-Binomial X
// (Double-Woodbury and Low-Rank Optimized)

functions {
    /**
     * MVN log density for covariance: Sigma = sigma_sq * I + Q Q'
     * using only sufficient statistics.
     *
     * Uses matrix determinant lemma and Woodbury identity:
     *   log|Sigma| = n*log(sigma_sq) + log|I + Q'Q/sigma_sq|
     *   r'Sigma^{-1}r = r'r/sigma_sq - (Q'r)'M^{-1}(Q'r)/sigma_sq^2
     *   where M = I + Q'Q/sigma_sq
     *
     * Complexity: O(k^3) for Cholesky of k×k, with k=2.
     */
    real mvn_diag_plus_lowrank_lpdf(real r2,
                                    vector Qtr,
                                    matrix QQ,
                                    real sigma_sq,
                                    int n) {
        int k = rows(Qtr);

        matrix[k, k] M = QQ / sigma_sq;
        for (i in 1:k) {
            M[i, i] += 1.0;
        }

        matrix[k, k] L_M = cholesky_decompose(M);
        real log_det = n * log(sigma_sq) + 2.0 * sum(log(diagonal(L_M)));

        vector[k] sol = mdivide_left_tri_low(L_M, Qtr);
        real quad = (r2 - dot_self(sol) / sigma_sq) / sigma_sq;

        return -0.5 * (n * log(2.0 * pi()) + log_det + quad);
    }

    /**
     * Compute log PMF of Beta-Binomial using recurrence relation
     *
     * BetaBinomial(k | n, alpha, beta) = C(n,k) * B(k+alpha, n-k+beta) / B(alpha, beta)
     * Recurrence: P(k)/P(k-1) = (n-k+1)(k-1+alpha) / [k(n-k+beta)]
     */
    vector beta_binomial_log_pmf_vec(int n, real alpha, real beta_param) {
        vector[n + 1] log_pmf;

        // Base case: k = 0
        log_pmf[1] = lgamma(n + beta_param) + lgamma(alpha + beta_param)
                     - lgamma(n + alpha + beta_param) - lgamma(beta_param);

        // Recurrence for k = 1, 2, ..., n
        for (k in 1:n) {
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

    // Number of trials for the Beta-Binomial
    array[G] int<lower=1> n_trials;

    array[N] int<lower=1, upper=G> id;      // Subject index per observation

    // Subject start position and length (for ragged array access)
    array[G] int<lower=1> pos;              // Start index for subject g
    array[G] int<lower=1> len;              // Number of obs for subject g
}

transformed data {
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

    // PRECOMPUTATION OPTIMIZATION: Compute constant time sums once
    vector[G] sumt;
    vector[G] sumt2;
    for (g in 1:G) {
        sumt[g] = 0.0;
        sumt2[g] = 0.0;
        int start_idx = pos[g];
        int n_g = len[g];
        for (i in 1:n_g) {
            real ti = t[start_idx + i - 1];
            sumt[g] += ti;
            sumt2[g] += ti * ti;
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
}

transformed parameters {
    // Cholesky factor of RE covariance: L_Sigma = diag(sigma_re) * L_corr
    matrix[2, 2] L_Sigma = diag_pre_multiply(sigma_re, L_re);
    real sigma_sq = square(sigma_main);

    // Imputation model: logit(mu) = linear predictor
    vector[G] eta_imputation = alpha_imputation + Z * gamma;
    vector[G] mu_imputation = inv_logit(eta_imputation);

    // Beta distribution parameters for each subject
    vector<lower=0>[G] alpha_beta = mu_imputation * phi_imputation;
    vector<lower=0>[G] beta_beta = (1 - mu_imputation) * phi_imputation;
    vector[N] xb = X * beta;

    // L_Sigma = [[a, 0], [b, c]]
    real re_a = L_Sigma[1, 1];
    real re_b = L_Sigma[2, 1];
    real re_c = L_Sigma[2, 2];
}

model {
    // Priors
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
    phi_imputation ~ gamma(2, 0.1);

    // =========================================================================
    // Likelihood for Subjects with OBSERVED x (O(n) Woodbury)
    // =========================================================================
    // Vectorize imputation likelihood for observed x
    x_obs ~ beta_binomial(m_trials, alpha_beta[index_obs], beta_beta[index_obs]);

    for (k in 1:G_obs) {
        int g = index_obs[k];
        int x_g = x_obs[k];
        int n_g = len[g];
        int start_idx = pos[g];

        real sum_r = 0;
        real sum_tr = 0;
        real sum_r2 = 0;

        for (i in 1:n_g) {
            int obs = start_idx + i - 1;
            real ti = t[obs];
            real mu_base_i = alpha_main + beta_t * ti + xb[obs];
            real ri = y[obs] - mu_base_i;

            sum_r += ri;
            sum_tr += ti * ri;
            sum_r2 += ri * ri;
        }

        real dd = n_g * square(beta_x)
                  + 2.0 * beta_x * beta_x_t_interaction * sumt[g]
                  + square(beta_x_t_interaction) * sumt2[g];
        real d_rbase = beta_x * sum_r + beta_x_t_interaction * sum_tr;
        real r2 = sum_r2 - 2.0 * x_g * d_rbase + square(x_g) * dd;

        matrix[2, 2] QQ_re;
        QQ_re[1, 1] = n_g * square(re_a) + 2.0 * re_a * re_b * sumt[g] + square(re_b) * sumt2[g];
        QQ_re[1, 2] = re_c * re_a * sumt[g] + re_c * re_b * sumt2[g];
        QQ_re[2, 1] = QQ_re[1, 2];
        QQ_re[2, 2] = square(re_c) * sumt2[g];

        vector[2] Qd;
        Qd[1] = n_g * re_a * beta_x
                + (re_a * beta_x_t_interaction + re_b * beta_x) * sumt[g]
                + re_b * beta_x_t_interaction * sumt2[g];
        Qd[2] = re_c * beta_x * sumt[g] + re_c * beta_x_t_interaction * sumt2[g];

        vector[2] Qtr_base;
        Qtr_base[1] = re_a * sum_r + re_b * sum_tr;
        Qtr_base[2] = re_c * sum_tr;

        vector[2] Qtr = Qtr_base - x_g * Qd;

        target += mvn_diag_plus_lowrank_lpdf(r2 | Qtr, QQ_re, sigma_sq, n_g);
    }

    // =========================================================================
    // Likelihood for Subjects with MISSING x (Fully Marginalized O(1) Woodbury Loop)
    // =========================================================================
    for (k in 1:G_mis) {
        int g = index_mis[k];
        int n_g = len[g];
        int start_idx = pos[g];

        real sum_r = 0;
        real sum_tr = 0;
        real sum_r2 = 0;

        for (i in 1:n_g) {
            int obs = start_idx + i - 1;
            real ti = t[obs];
            real mu_base_i = alpha_main + beta_t * ti + xb[obs];
            real ri = y[obs] - mu_base_i;

            sum_r += ri;
            sum_tr += ti * ri;
            sum_r2 += ri * ri;
        }

        real dd = n_g * square(beta_x)
                  + 2.0 * beta_x * beta_x_t_interaction * sumt[g]
                  + square(beta_x_t_interaction) * sumt2[g];
        real d_rbase = beta_x * sum_r + beta_x_t_interaction * sum_tr;

        matrix[2, 2] QQ_re;
        QQ_re[1, 1] = n_g * square(re_a) + 2.0 * re_a * re_b * sumt[g] + square(re_b) * sumt2[g];
        QQ_re[1, 2] = re_c * re_a * sumt[g] + re_c * re_b * sumt2[g];
        QQ_re[2, 1] = QQ_re[1, 2];
        QQ_re[2, 2] = square(re_c) * sumt2[g];

        vector[2] Qd;
        Qd[1] = n_g * re_a * beta_x
                + (re_a * beta_x_t_interaction + re_b * beta_x) * sumt[g]
                + re_b * beta_x_t_interaction * sumt2[g];
        Qd[2] = re_c * beta_x * sumt[g] + re_c * beta_x_t_interaction * sumt2[g];

        vector[2] Qtr_base;
        Qtr_base[1] = re_a * sum_r + re_b * sum_tr;
        Qtr_base[2] = re_c * sum_tr;

        // Compute all Beta-Binomial log PMFs using efficient recurrence
        vector[m_trials + 1] log_pmf_bb = beta_binomial_log_pmf_vec(
            m_trials, alpha_beta[g], beta_beta[g]);

        // Precompute Woodbury linear algebra outside the inner loop
        matrix[2, 2] M = QQ_re / sigma_sq;
        M[1, 1] += 1.0;
        M[2, 2] += 1.0;
        matrix[2, 2] L_M = cholesky_decompose(M);

        real log_det = n_g * log(sigma_sq) + 2.0 * sum(log(diagonal(L_M)));
        real mvn_const = -0.5 * (n_g * log(2.0 * pi()) + log_det);

        vector[2] a_w = mdivide_left_tri_low(L_M, Qtr_base);
        vector[2] b_w = mdivide_left_tri_low(L_M, Qd);

        // Collapse Woodbury quadratic forms into a single j-polynomial
        // so the inner loop matches main's per-iter autodiff cost (~5 ops vs ~12).
        real inv_s = 1.0 / sigma_sq;
        real inv_s2 = inv_s * inv_s;
        real K1 = sum_r2 * inv_s - dot_self(a_w) * inv_s2;
        real K2 = d_rbase * inv_s - dot_product(a_w, b_w) * inv_s2;
        real K3 = dd * inv_s - dot_self(b_w) * inv_s2;

        vector[m_trials + 1] log_probs;
        for (j in 0:m_trials) {
            real Q = K1 - 2.0 * j * K2 + square(j) * K3;
            real lp_y_xj = mvn_const - 0.5 * Q;
            log_probs[j + 1] = lp_y_xj + log_pmf_bb[j + 1];
        }

        target += log_sum_exp(log_probs);
    }
}

generated quantities {
    corr_matrix[2] corr_rand_effects = multiply_lower_tri_self_transpose(L_re);
    cov_matrix[2] cov_rand_effects = quad_form_diag(corr_rand_effects, sigma_re);

    vector[G] log_lik;

    for (k in 1:G_obs) {
        int g = index_obs[k];
        int x_g = x_obs[k];
        int n_g = len[g];
        int start_idx = pos[g];

        real sum_r = 0;
        real sum_tr = 0;
        real sum_r2 = 0;

        for (i in 1:n_g) {
            int obs = start_idx + i - 1;
            real ti = t[obs];
            real ri = y[obs] - (alpha_main + beta_t * ti + xb[obs]);
            sum_r += ri;
            sum_tr += ti * ri;
            sum_r2 += ri * ri;
        }

        real dd = n_g * square(beta_x)
                  + 2.0 * beta_x * beta_x_t_interaction * sumt[g]
                  + square(beta_x_t_interaction) * sumt2[g];
        real d_rbase = beta_x * sum_r + beta_x_t_interaction * sum_tr;
        real r2 = sum_r2 - 2.0 * x_g * d_rbase + square(x_g) * dd;

        matrix[2, 2] QQ_re;
        QQ_re[1, 1] = n_g * square(re_a) + 2.0 * re_a * re_b * sumt[g] + square(re_b) * sumt2[g];
        QQ_re[1, 2] = re_c * re_a * sumt[g] + re_c * re_b * sumt2[g];
        QQ_re[2, 1] = QQ_re[1, 2];
        QQ_re[2, 2] = square(re_c) * sumt2[g];

        vector[2] Qd;
        Qd[1] = n_g * re_a * beta_x
                + (re_a * beta_x_t_interaction + re_b * beta_x) * sumt[g]
                + re_b * beta_x_t_interaction * sumt2[g];
        Qd[2] = re_c * beta_x * sumt[g] + re_c * beta_x_t_interaction * sumt2[g];

        vector[2] Qtr_base;
        Qtr_base[1] = re_a * sum_r + re_b * sum_tr;
        Qtr_base[2] = re_c * sum_tr;

        vector[2] Qtr = Qtr_base - x_g * Qd;

        log_lik[g] = mvn_diag_plus_lowrank_lpdf(r2 | Qtr, QQ_re, sigma_sq, n_g)
                     + beta_binomial_lpmf(x_g | m_trials, alpha_beta[g], beta_beta[g]);
    }

    for (k in 1:G_mis) {
        int g = index_mis[k];
        int n_g = len[g];
        int start_idx = pos[g];

        real sum_r = 0;
        real sum_tr = 0;
        real sum_r2 = 0;

        for (i in 1:n_g) {
            int obs = start_idx + i - 1;
            real ti = t[obs];
            real ri = y[obs] - (alpha_main + beta_t * ti + xb[obs]);
            sum_r += ri;
            sum_tr += ti * ri;
            sum_r2 += ri * ri;
        }

        real dd = n_g * square(beta_x)
                  + 2.0 * beta_x * beta_x_t_interaction * sumt[g]
                  + square(beta_x_t_interaction) * sumt2[g];
        real d_rbase = beta_x * sum_r + beta_x_t_interaction * sum_tr;

        matrix[2, 2] QQ_re;
        QQ_re[1, 1] = n_g * square(re_a) + 2.0 * re_a * re_b * sumt[g] + square(re_b) * sumt2[g];
        QQ_re[1, 2] = re_c * re_a * sumt[g] + re_c * re_b * sumt2[g];
        QQ_re[2, 1] = QQ_re[1, 2];
        QQ_re[2, 2] = square(re_c) * sumt2[g];

        vector[2] Qd;
        Qd[1] = n_g * re_a * beta_x
                + (re_a * beta_x_t_interaction + re_b * beta_x) * sumt[g]
                + re_b * beta_x_t_interaction * sumt2[g];
        Qd[2] = re_c * beta_x * sumt[g] + re_c * beta_x_t_interaction * sumt2[g];

        vector[2] Qtr_base;
        Qtr_base[1] = re_a * sum_r + re_b * sum_tr;
        Qtr_base[2] = re_c * sum_tr;

        vector[m_trials + 1] log_pmf_bb = beta_binomial_log_pmf_vec(
            m_trials, alpha_beta[g], beta_beta[g]);

        vector[m_trials + 1] log_probs;
        for (j in 0:m_trials) {
            real r2_j = sum_r2 - 2.0 * j * d_rbase + square(j) * dd;
            vector[2] Qtr_j = Qtr_base - j * Qd;

            real lp_y_xj = mvn_diag_plus_lowrank_lpdf(r2_j | Qtr_j, QQ_re, sigma_sq, n_g);
            log_probs[j + 1] = lp_y_xj + log_pmf_bb[j + 1];
        }

        log_lik[g] = log_sum_exp(log_probs);
    }
}
