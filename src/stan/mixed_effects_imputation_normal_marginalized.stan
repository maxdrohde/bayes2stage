// Mixed Effects Model with Fully Marginalized Random Effects AND Missing X
// (Double-Woodbury Optimized with Precomputed Sums)

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
     * Complexity: O(k^3) for Cholesky of k×k, with k=2 or 3.
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
}

data {
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

    array[G_obs] int<lower=1, upper=G> index_obs;
    array[G_mis] int<lower=1, upper=G> index_mis;
    vector[G_obs] x_obs;

    array[G] int<lower=1> pos;
    array[G] int<lower=1> len;
}

transformed data {
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
    vector<lower=0>[2] sigma_re;
    cholesky_factor_corr[2] L_re;

    real alpha_main;
    real beta_t;
    real beta_x;
    real beta_x_t_interaction;
    vector[P] beta;
    real<lower=1e-8> sigma_main;

    real alpha_imputation;
    vector[S] gamma;
    real<lower=0> sigma_imputation;
}

transformed parameters {
    matrix[2, 2] L_Sigma = diag_pre_multiply(sigma_re, L_re);
    real sigma_sq = square(sigma_main);
    vector[N] xb = X * beta;
    vector[G] x_imputation_mean = alpha_imputation + Z * gamma;

    // L_Sigma = [[a, 0], [b, c]]
    real re_a = L_Sigma[1, 1];
    real re_b = L_Sigma[2, 1];
    real re_c = L_Sigma[2, 2];
}

model {
    // Priors
    sigma_re ~ exponential(0.1);
    L_re ~ lkj_corr_cholesky(2);
    sigma_main ~ exponential(0.1);
    alpha_main ~ normal(0, 100);
    beta_t ~ normal(0, 100);
    beta_x ~ normal(0, 100);
    beta_x_t_interaction ~ normal(0, 100);
    beta ~ normal(0, 100);
    sigma_imputation ~ exponential(0.1);
    alpha_imputation ~ normal(0, 100);
    gamma ~ normal(0, 100);

    // OBSERVED x: diagonal + rank-2
    for (k_subj in 1:G_obs) {
        int g = index_obs[k_subj];
        int n_g = len[g];
        int start_idx = pos[g];
        real x_g = x_obs[k_subj];

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
        x_g ~ normal(x_imputation_mean[g], sigma_imputation);
    }

    // MISSING x: diagonal + rank-3
    for (k_subj in 1:G_mis) {
        int g = index_mis[k_subj];
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

        real mu_x = x_imputation_mean[g];
        real sx = sigma_imputation;

        real dd = n_g * square(beta_x)
                  + 2.0 * beta_x * beta_x_t_interaction * sumt[g]
                  + square(beta_x_t_interaction) * sumt2[g];
        real d_rbase = beta_x * sum_r + beta_x_t_interaction * sum_tr;
        real r2 = sum_r2 - 2.0 * mu_x * d_rbase + square(mu_x) * dd;

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

        vector[2] Qtr_re = Qtr_base - mu_x * Qd;
        real d_r = d_rbase - mu_x * dd;

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

        target += mvn_diag_plus_lowrank_lpdf(r2 | Qtr_full, QQ_full, sigma_sq, n_g);
    }
}
