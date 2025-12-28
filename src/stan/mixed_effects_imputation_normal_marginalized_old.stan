// OLD VERSION: Mixed Effects Model with Marginalized RE and Missing X
// Uses O(n³) chol_rank1_update instead of O(n²) Woodbury identity
// This file is for benchmarking comparison only

functions {
    matrix compute_marginal_cov_chol(vector t_g, matrix L_Sigma, real sigma) {
        int n = rows(t_g);
        matrix[n, 2] Z_re;
        matrix[n, n] V;

        Z_re[, 1] = rep_vector(1.0, n);
        Z_re[, 2] = t_g;

        V = tcrossprod(Z_re * L_Sigma);

        for (i in 1:n) {
            V[i, i] += square(sigma);
        }

        return cholesky_decompose(V);
    }

    // OLD: O(n³) Cholesky update
    matrix chol_rank1_update(matrix L_V, vector d, real sigma_x) {
        int n = rows(L_V);
        matrix[n, n] V = tcrossprod(L_V);

        V += square(sigma_x) * tcrossprod(to_matrix(d, n, 1));

        return cholesky_decompose(V);
    }
}

data {
    int<lower=0> N;
    int<lower=1> G;
    int<lower=0, upper=G> G_obs;
    int<lower=0, upper=G> G_mis;

    int<lower=0> P;
    vector[N] t;
    matrix[N, P] X;

    int<lower=0> S;
    matrix[G, S] Z;

    vector[N] y;

    array[G_obs] int<lower=1, upper=G> index_obs;
    array[G_mis] int<lower=1, upper=G> index_mis;
    vector[G_obs] x_obs;

    array[N] int<lower=1, upper=G> id;

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
    real<lower=0> sigma_main;

    real alpha_imputation;
    vector[S] gamma;
    real<lower=0> sigma_imputation;
}

transformed parameters {
    matrix[2, 2] L_Sigma = diag_pre_multiply(sigma_re, L_re);
    vector[G] x_imputation_mean = alpha_imputation + Z * gamma;
}

model {
    sigma_re ~ exponential(0.1);
    L_re ~ lkj_corr_cholesky(2);

    sigma_main ~ exponential(0.1);
    alpha_main ~ normal(0, 100);
    beta ~ normal(0, 100);
    beta_t ~ normal(0, 100);
    beta_x ~ normal(0, 100);
    beta_x_t_interaction ~ normal(0, 100);

    sigma_imputation ~ exponential(0.1);
    alpha_imputation ~ normal(0, 100);
    gamma ~ normal(0, 100);

    // Likelihood for Subjects with OBSERVED x
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

        y_g ~ multi_normal_cholesky(mu_g, L_V);
        x_g ~ normal(x_imputation_mean[g], sigma_imputation);
    }

    // Likelihood for Subjects with MISSING x (OLD: O(n³) Cholesky update)
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

        // OLD: O(n³) rank-1 Cholesky update
        matrix[n_g, n_g] L_V_aug = chol_rank1_update(L_V, d_g, sigma_imputation);

        y_g ~ multi_normal_cholesky(mu_g, L_V_aug);
    }
}

// generated quantities {
//     corr_matrix[2] corr_rand_effects = multiply_lower_tri_self_transpose(L_re);
//     cov_matrix[2] cov_rand_effects = quad_form_diag(corr_rand_effects, sigma_re);
//
//     vector[G] log_lik;
//
//     for (k in 1:G_obs) {
//         int g = index_obs[k];
//         real x_g = x_obs[k];
//         int n_g = len[g];
//         int start_idx = pos[g];
//
//         vector[n_g] y_g = segment(y, start_idx, n_g);
//         vector[n_g] t_g = segment(t, start_idx, n_g);
//
//         vector[n_g] c_g;
//         for (i in 1:n_g) {
//             int obs_idx = start_idx + i - 1;
//             c_g[i] = alpha_main + beta_t * t_g[i] + dot_product(X[obs_idx], beta);
//         }
//         vector[n_g] d_g = rep_vector(beta_x, n_g) + beta_x_t_interaction * t_g;
//         vector[n_g] mu_g = c_g + x_g * d_g;
//
//         matrix[n_g, n_g] L_V = compute_marginal_cov_chol(t_g, L_Sigma, sigma_main);
//
//         real lp_y = multi_normal_cholesky_lpdf(y_g | mu_g, L_V);
//         real lp_x = normal_lpdf(x_g | x_imputation_mean[g], sigma_imputation);
//
//         log_lik[g] = lp_y + lp_x;
//     }
//
//     for (k in 1:G_mis) {
//         int g = index_mis[k];
//         int n_g = len[g];
//         int start_idx = pos[g];
//
//         vector[n_g] y_g = segment(y, start_idx, n_g);
//         vector[n_g] t_g = segment(t, start_idx, n_g);
//
//         vector[n_g] c_g;
//         for (i in 1:n_g) {
//             int obs_idx = start_idx + i - 1;
//             c_g[i] = alpha_main + beta_t * t_g[i] + dot_product(X[obs_idx], beta);
//         }
//         vector[n_g] d_g = rep_vector(beta_x, n_g) + beta_x_t_interaction * t_g;
//
//         real mu_x = x_imputation_mean[g];
//         vector[n_g] mu_g = c_g + mu_x * d_g;
//
//         matrix[n_g, n_g] L_V = compute_marginal_cov_chol(t_g, L_Sigma, sigma_main);
//         matrix[n_g, n_g] L_V_aug = chol_rank1_update(L_V, d_g, sigma_imputation);
//
//         log_lik[g] = multi_normal_cholesky_lpdf(y_g | mu_g, L_V_aug);
//     }
// }
