// Mixed Effects Model with Fully Marginalized Random Effects AND Missing X
// (Readable version — mathematically equivalent to mixed_effects_imputation_normal_marginalized.stan)
//
// MATHEMATICAL SETUP
// ------------------
// Outcome model (linear mixed effects):
//   y_ig = alpha + beta_x * x_i + beta_t * t_ig + beta_xt * x_i * t_ig + X_ig * beta
//          + b_0i + b_1i * t_ig + eps_ig
//
//   where (b_0i, b_1i) ~ MVN(0, Sigma_re) are correlated random intercepts/slopes
//   and eps_ig ~ N(0, sigma^2)
//
// Imputation model (for missing x):
//   x_i ~ N(alpha_imp + Z_i * gamma, sigma_imp^2)
//
// MARGINALIZATION STRATEGY
// -----------------------
// Integrating out random effects analytically:
//   y_i | x_i ~ MVN(mu_i, V_i)
//   V_i = sigma^2 * I + Z_i * Sigma_re * Z_i'   (diagonal + rank-2)
//
// For subjects with missing x, additionally integrate out x_i:
//   y_i ~ MVN(mu_i(E[x_i]), V_i + sigma_x^2 * d_i * d_i')   (diagonal + rank-3)
//   where d_ig = beta_x + beta_xt * t_ig is the coefficient of x in observation ig.
//
// WOODBURY IDENTITY
// -----------------
// V_i has structure: sigma^2 * I + Q * Q'  where Q = Z_i * L_Sigma (the Cholesky root).
// The Woodbury identity evaluates log|V| and r' V^{-1} r using only k x k
// operations (k = 2 for observed x, k = 3 for missing x), avoiding O(n^3) cost.
//
// Result: ZERO latent variables. Complexity O(n) + O(k^3) per subject.

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
     * Complexity: O(k^3) for Cholesky of k x k, with k=2 or 3.
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
     * Compute sufficient statistics for one subject's observations.
     *
     * Given time points t and base residuals r = y - (alpha + beta_t * t + X * beta),
     * computes the 5 scalar summaries needed for the Woodbury identity:
     *   sum_t   = sum(t_ig)
     *   sum_t2  = sum(t_ig^2)
     *   sum_r   = sum(r_ig)        where r_ig are base residuals (excluding x terms)
     *   sum_tr  = sum(t_ig * r_ig)
     *   sum_r2  = sum(r_ig^2)
     *
     * Returns: vector[5] = {sum_t, sum_t2, sum_r, sum_tr, sum_r2}
     */
    vector compute_sufficient_stats(vector t_subj, vector y_subj,
                                    real alpha_main, real beta_t,
                                    vector xb_subj, int n_g) {
        real sum_t = 0;
        real sum_t2 = 0;
        real sum_r = 0;
        real sum_tr = 0;
        real sum_r2 = 0;

        for (i in 1:n_g) {
            real ti = t_subj[i];
            real ri = y_subj[i] - (alpha_main + beta_t * ti + xb_subj[i]);

            sum_t += ti;
            sum_t2 += ti * ti;
            sum_r += ri;
            sum_tr += ti * ri;
            sum_r2 += ri * ri;
        }

        return [sum_t, sum_t2, sum_r, sum_tr, sum_r2]';
    }

    /**
     * Build the 2x2 random-effects Gram matrix Q'Q.
     *
     * The random-effects design for subject i is Z_i = [1, t_i1; 1, t_i2; ...] (n_g x 2).
     * Multiplying by L_Sigma (Cholesky of Sigma_re) gives Q = Z_i * L_Sigma (n_g x 2).
     * The Gram matrix Q'Q is 2x2 and expressible via sufficient statistics:
     *
     *   Q'Q = [ L11^2*n + 2*L11*L21*sum_t + L21^2*sum_t2,   L22*(L11*sum_t + L21*sum_t2) ]
     *         [              (symmetric)                  ,   L22^2*sum_t2                  ]
     */
    matrix build_re_gram_matrix(real L11, real L21, real L22,
                                int n_g, real sum_t, real sum_t2) {
        matrix[2, 2] QQ;
        QQ[1, 1] = n_g * square(L11) + 2.0 * L11 * L21 * sum_t + square(L21) * sum_t2;
        QQ[1, 2] = L22 * L11 * sum_t + L22 * L21 * sum_t2;
        QQ[2, 1] = QQ[1, 2];
        QQ[2, 2] = square(L22) * sum_t2;
        return QQ;
    }

    /**
     * Compute Q' * r_base (projection of base residuals onto random-effects space).
     *
     * Q_resid = (Z_i * L_Sigma)' * r_base
     *         = [ L11 * sum_r + L21 * sum_tr,  L22 * sum_tr ]'
     */
    vector compute_Q_resid(real L11, real L21, real L22,
                           real sum_r, real sum_tr) {
        return [L11 * sum_r + L21 * sum_tr, L22 * sum_tr]';
    }

    /**
     * Compute how x enters the residuals via the coefficient vector d.
     *
     * The mean function includes x_i via: beta_x * x_i + beta_xt * x_i * t_ig.
     * Define d_ig = beta_x + beta_xt * t_ig (the coefficient of x in observation ig).
     * Then:
     *   x_quad_coeff  = d'd = sum(d_ig^2)     -- quadratic coefficient of x in ||r||^2
     *   x_resid_cross = d' r_base             -- cross term between x and base residuals
     *   Q_x_coeff     = Q'd                   -- projection of d onto random-effects space
     *
     * Returns: vector[4] = {x_quad_coeff, x_resid_cross, Q_x_coeff[1], Q_x_coeff[2]}
     */
    vector compute_x_coefficients(real beta_x, real beta_xt,
                                  real L11, real L21, real L22,
                                  int n_g, real sum_t, real sum_t2,
                                  real sum_r, real sum_tr) {
        // d'd = ||d||^2 where d_ig = beta_x + beta_xt * t_ig
        real x_quad_coeff = n_g * square(beta_x)
                            + 2.0 * beta_x * beta_xt * sum_t
                            + square(beta_xt) * sum_t2;

        // d' * r_base = sum(d_ig * r_ig)
        real x_resid_cross = beta_x * sum_r + beta_xt * sum_tr;

        // Q' * d = (Z * L_Sigma)' * d
        real Q_x_coeff_1 = n_g * L11 * beta_x
                           + (L11 * beta_xt + L21 * beta_x) * sum_t
                           + L21 * beta_xt * sum_t2;
        real Q_x_coeff_2 = L22 * beta_x * sum_t + L22 * beta_xt * sum_t2;

        return [x_quad_coeff, x_resid_cross, Q_x_coeff_1, Q_x_coeff_2]';
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

    // Cholesky factor of random-effects covariance:
    //   L_Sigma = [[L11,   0],
    //              [L21, L22]]
    real L11 = L_Sigma[1, 1];
    real L21 = L_Sigma[2, 1];
    real L22 = L_Sigma[2, 2];
}

model {
    // =========================================================================
    // PRIORS
    // =========================================================================
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

    // =========================================================================
    // OBSERVED x: marginal likelihood via rank-2 Woodbury
    // =========================================================================
    // For each subject with observed x_i:
    //   y_i | x_i ~ MVN(mu_i, sigma^2 I + Q Q')   where Q = Z_i * L_Sigma
    //
    // Steps:
    //   1. Compute sufficient statistics from base residuals (excluding x terms)
    //   2. Compute x-dependent coefficients (how x enters residuals)
    //   3. Form ||r||^2, Q'Q (2x2), and Q'r (2x1) for the full residuals
    //   4. Evaluate MVN density via rank-2 Woodbury

    for (k_subj in 1:G_obs) {
        int g = index_obs[k_subj];
        int n_g = len[g];
        int start_idx = pos[g];
        real x_g = x_obs[k_subj];

        // Step 1: sufficient statistics from base residuals
        vector[5] ss = compute_sufficient_stats(
            t[start_idx:(start_idx + n_g - 1)],
            y[start_idx:(start_idx + n_g - 1)],
            alpha_main, beta_t,
            xb[start_idx:(start_idx + n_g - 1)],
            n_g
        );
        real sum_t  = ss[1];
        real sum_t2 = ss[2];
        real sum_r  = ss[3];
        real sum_tr = ss[4];
        real sum_r2 = ss[5];

        // Step 2: how x enters the residuals
        //   d_ig = beta_x + beta_xt * t_ig  (coefficient of x_i in observation ig)
        vector[4] xc = compute_x_coefficients(
            beta_x, beta_x_t_interaction,
            L11, L21, L22,
            n_g, sum_t, sum_t2, sum_r, sum_tr
        );
        real x_quad_coeff       = xc[1];   // d'd
        real x_resid_cross      = xc[2];   // d' r_base
        vector[2] Q_x_coeff     = xc[3:4]; // Q'd

        // Step 3: full residual norm and projections
        //   Full residual: r = r_base - x_g * d
        //   ||r||^2 = ||r_base||^2 - 2*x_g*(d'r_base) + x_g^2*(d'd)
        real r2 = sum_r2 - 2.0 * x_g * x_resid_cross + square(x_g) * x_quad_coeff;

        matrix[2, 2] QQ = build_re_gram_matrix(L11, L21, L22, n_g, sum_t, sum_t2);
        vector[2] Q_resid = compute_Q_resid(L11, L21, L22, sum_r, sum_tr);
        //   Q'r = Q'r_base - x_g * Q'd
        vector[2] Qtr = Q_resid - x_g * Q_x_coeff;

        // Step 4: marginal log-likelihood via Woodbury (k=2)
        target += mvn_diag_plus_lowrank_lpdf(r2 | Qtr, QQ, sigma_sq, n_g);

        // Imputation model likelihood for observed x
        x_g ~ normal(x_imputation_mean[g], sigma_imputation);
    }

    // =========================================================================
    // MISSING x: marginal likelihood via rank-3 Woodbury
    // =========================================================================
    // For each subject with missing x_i, additionally integrate out x_i:
    //   y_i ~ MVN(mu_i(E[x_i]), sigma^2 I + Q_aug Q_aug')
    //
    // The augmented Q matrix has 3 columns:
    //   Columns 1-2: Z_i * L_Sigma     (random effects, same as observed case)
    //   Column 3:    sigma_x * d_i      (x uncertainty; d_ig = beta_x + beta_xt * t_ig)
    //
    // This adds one rank to the covariance, capturing the additional variance
    // from not knowing x_i.

    for (k_subj in 1:G_mis) {
        int g = index_mis[k_subj];
        int n_g = len[g];
        int start_idx = pos[g];

        real mu_x = x_imputation_mean[g];
        real sigma_x = sigma_imputation;

        // Step 1: sufficient statistics (same as observed-x case)
        vector[5] ss = compute_sufficient_stats(
            t[start_idx:(start_idx + n_g - 1)],
            y[start_idx:(start_idx + n_g - 1)],
            alpha_main, beta_t,
            xb[start_idx:(start_idx + n_g - 1)],
            n_g
        );
        real sum_t  = ss[1];
        real sum_t2 = ss[2];
        real sum_r  = ss[3];
        real sum_tr = ss[4];
        real sum_r2 = ss[5];

        // Step 2: x-dependent coefficients (same computation)
        vector[4] xc = compute_x_coefficients(
            beta_x, beta_x_t_interaction,
            L11, L21, L22,
            n_g, sum_t, sum_t2, sum_r, sum_tr
        );
        real x_quad_coeff       = xc[1];   // d'd
        real x_resid_cross      = xc[2];   // d' r_base
        vector[2] Q_x_coeff     = xc[3:4]; // Q'd

        // Step 3a: residual norm evaluated at x = mu_x (the imputation mean)
        real r2 = sum_r2 - 2.0 * mu_x * x_resid_cross + square(mu_x) * x_quad_coeff;

        // Step 3b: rank-2 components (same as observed, but evaluated at mu_x)
        matrix[2, 2] QQ_re = build_re_gram_matrix(L11, L21, L22, n_g, sum_t, sum_t2);
        vector[2] Q_resid = compute_Q_resid(L11, L21, L22, sum_r, sum_tr);
        vector[2] Qtr_re = Q_resid - mu_x * Q_x_coeff;

        // Cross term between d and residual at mu_x: d'(r_base - mu_x * d)
        real x_adjusted_cross = x_resid_cross - mu_x * x_quad_coeff;

        // Step 3c: augment to rank-3 by adding x-uncertainty column
        //   Q_aug = [Q_re | sigma_x * d]   (n_g x 3 matrix)
        //
        //   Q_aug' Q_aug = [ QQ_re          , sigma_x * Q'd     ]
        //                  [ sigma_x * (Q'd)', sigma_x^2 * d'd   ]
        matrix[3, 3] QQ_aug;
        QQ_aug[1:2, 1:2] = QQ_re;
        QQ_aug[1, 3] = sigma_x * Q_x_coeff[1];
        QQ_aug[2, 3] = sigma_x * Q_x_coeff[2];
        QQ_aug[3, 1] = sigma_x * Q_x_coeff[1];
        QQ_aug[3, 2] = sigma_x * Q_x_coeff[2];
        QQ_aug[3, 3] = square(sigma_x) * x_quad_coeff;

        //   Q_aug' r = [ Q_re' r                       ]
        //              [ sigma_x * d'(r_base - mu_x*d)  ]
        vector[3] Qtr_aug;
        Qtr_aug[1:2] = Qtr_re;
        Qtr_aug[3] = sigma_x * x_adjusted_cross;

        // Step 4: marginal log-likelihood via Woodbury (k=3)
        target += mvn_diag_plus_lowrank_lpdf(r2 | Qtr_aug, QQ_aug, sigma_sq, n_g);
    }
}

