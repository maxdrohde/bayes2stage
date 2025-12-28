# Compare double-Woodbury model to standard mixture model
# Timing and parameter recovery comparison

library(cmdstanr)
library(dplyr)
library(tidyr)
library(tibble)

set.seed(777)

# ==============================================================================
# Simulate data
# ==============================================================================

# Settings
G <- 200          # Number of subjects
n_per_subject <- 10  # Observations per subject (larger = more speedup expected)
missing_prop <- 0.3  # Proportion with missing x
K <- 3            # Mixture components

# True parameters
TRUE_PARAMS <- list(
    alpha_main = 5,
    beta_t = 0.5,
    beta_x = 2,
    beta_x_t_interaction = 0.3,
    sigma_main = 1.5,
    sigma_re_intercept = 2,
    sigma_re_slope = 0.5,
    rho_re = 0.3,
    alpha_imputation = 0,
    sigma_component = c(0.8, 1.0, 1.2),
    mu_component = c(-2, 0, 2),
    theta = c(0.3, 0.4, 0.3)
)

# Generate subject-level data
subjects <- tibble(
    subject_id = 1:G,
    has_missing_x = sample(c(TRUE, FALSE), G, replace = TRUE,
                           prob = c(missing_prop, 1 - missing_prop))
)

# Generate x from mixture
generate_x_mixture <- function(n, theta, mu, sigma) {
    components <- sample(1:length(theta), n, replace = TRUE, prob = theta)
    rnorm(n, mean = mu[components], sd = sigma[components])
}

subjects <- subjects |>
    mutate(
        x_true = generate_x_mixture(G, TRUE_PARAMS$theta,
                                    TRUE_PARAMS$mu_component,
                                    TRUE_PARAMS$sigma_component),
        x_obs = if_else(has_missing_x, NA_real_, x_true)
    )

# Generate random effects
Sigma_re <- matrix(c(
    TRUE_PARAMS$sigma_re_intercept^2,
    TRUE_PARAMS$rho_re * TRUE_PARAMS$sigma_re_intercept * TRUE_PARAMS$sigma_re_slope,
    TRUE_PARAMS$rho_re * TRUE_PARAMS$sigma_re_intercept * TRUE_PARAMS$sigma_re_slope,
    TRUE_PARAMS$sigma_re_slope^2
), nrow = 2)

L_re <- t(chol(Sigma_re))
re <- matrix(rnorm(G * 2), nrow = G) %*% t(L_re)
subjects$re_intercept <- re[, 1]
subjects$re_slope <- re[, 2]

# Generate observation-level data
obs_data <- subjects |>
    slice(rep(1:n(), each = n_per_subject)) |>
    group_by(subject_id) |>
    mutate(
        obs_id = row_number(),
        t = seq(0, 1, length.out = n_per_subject)[obs_id]
    ) |>
    ungroup() |>
    mutate(
        mu = TRUE_PARAMS$alpha_main +
             TRUE_PARAMS$beta_t * t +
             TRUE_PARAMS$beta_x * x_true +
             TRUE_PARAMS$beta_x_t_interaction * x_true * t +
             re_intercept + re_slope * t,
        y = rnorm(n(), mu, TRUE_PARAMS$sigma_main)
    )

# ==============================================================================
# Prepare Stan data
# ==============================================================================

# Subject indices
index_obs <- which(!subjects$has_missing_x)
index_mis <- which(subjects$has_missing_x)

# Ragged array info
pos <- seq(1, nrow(obs_data), by = n_per_subject)
len <- rep(n_per_subject, G)

stan_data <- list(
    N = nrow(obs_data),
    G = G,
    G_obs = length(index_obs),
    G_mis = length(index_mis),
    P = 0,
    t = obs_data$t,
    X = matrix(0, nrow = nrow(obs_data), ncol = 0),
    y = obs_data$y,
    S = 0,
    Z = matrix(0, nrow = G, ncol = 0),
    index_obs = index_obs,
    index_mis = index_mis,
    x_obs = subjects$x_obs[index_obs],
    id = obs_data$subject_id,
    pos = pos,
    len = len,
    K = K
)

# Fix zero-column matrices for Stan
stan_data$P <- 1
stan_data$X <- matrix(0, nrow = nrow(obs_data), ncol = 1)
stan_data$S <- 1
stan_data$Z <- matrix(0, nrow = G, ncol = 1)

cat("Data summary:\n")
cat("  Subjects:", G, "\n")
cat("  Obs per subject:", n_per_subject, "\n
")
cat("  Total obs:", nrow(obs_data), "\n")
cat("  Subjects with observed x:", length(index_obs), "\n")
cat("  Subjects with missing x:", length(index_mis), "\n")

# ==============================================================================
# Compile models
# ==============================================================================

cat("\nCompiling models...\n")

model_mixture <- cmdstan_model(
    "src/stan/mixed_effects_imputation_normal_marginalized_mixture.stan",
    quiet = TRUE
)

model_double_woodbury <- cmdstan_model(
    "src/stan/mixed_effects_imputation_normal_marginalized_double_woodbury.stan",
    quiet = TRUE
)

cat("Models compiled.\n")

# ==============================================================================
# Fit models
# ==============================================================================

n_chains <- 4
n_iter <- 1000
n_warmup <- 500

cat("\nFitting mixture model...\n")
time_mixture <- system.time({
    fit_mixture <- model_mixture$sample(
        data = stan_data,
        chains = n_chains,
        parallel_chains = n_chains,
        iter_warmup = n_warmup,
        iter_sampling = n_iter - n_warmup,
        refresh = 0,
        show_messages = FALSE
    )
})

cat("Fitting double-Woodbury model...\n")
time_double_woodbury <- system.time({
    fit_double_woodbury <- model_double_woodbury$sample(
        data = stan_data,
        chains = n_chains,
        parallel_chains = n_chains,
        iter_warmup = n_warmup,
        iter_sampling = n_iter - n_warmup,
        refresh = 0,
        show_messages = FALSE
    )
})

# ==============================================================================
# Compare timing
# ==============================================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("TIMING COMPARISON\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")

cat("\nMixture model:\n")
cat("  Elapsed time:", round(time_mixture["elapsed"], 2), "seconds\n")

cat("\nDouble-Woodbury model:\n")
cat("  Elapsed time:", round(time_double_woodbury["elapsed"], 2), "seconds\n")

speedup <- time_mixture["elapsed"] / time_double_woodbury["elapsed"]
cat("\nSpeedup:", round(speedup, 2), "x\n")

# ==============================================================================
# Compare parameter estimates
# ==============================================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("PARAMETER ESTIMATES\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")

params_to_compare <- c(
    "alpha_main", "beta_t", "beta_x", "beta_x_t_interaction",
    "sigma_main", "sigma_re[1]", "sigma_re[2]",
    "alpha_imputation", "theta[1]", "theta[2]", "theta[3]",
    "sigma_component[1]", "sigma_component[2]", "sigma_component[3]"
)

true_values <- c(
    TRUE_PARAMS$alpha_main,
    TRUE_PARAMS$beta_t,
    TRUE_PARAMS$beta_x,
    TRUE_PARAMS$beta_x_t_interaction,
    TRUE_PARAMS$sigma_main,
    TRUE_PARAMS$sigma_re_intercept,
    TRUE_PARAMS$sigma_re_slope,
    TRUE_PARAMS$alpha_imputation,
    TRUE_PARAMS$theta,
    TRUE_PARAMS$sigma_component
)

summary_mixture <- fit_mixture$summary(params_to_compare)
summary_dw <- fit_double_woodbury$summary(params_to_compare)

comparison <- tibble(
    parameter = params_to_compare,
    true_value = true_values,
    mixture_mean = summary_mixture$mean,
    mixture_sd = summary_mixture$sd,
    dw_mean = summary_dw$mean,
    dw_sd = summary_dw$sd,
    diff = abs(summary_mixture$mean - summary_dw$mean)
)

cat("\n")
cat(sprintf("%-25s %8s %12s %12s %10s\n",
            "Parameter", "True", "Mixture", "DoubleWood", "Diff"))
cat(rep("-", 70) |> paste(collapse = ""), "\n")

for (i in 1:nrow(comparison)) {
    cat(sprintf("%-25s %8.3f %7.3f (%.2f) %7.3f (%.2f) %10.4f\n",
                comparison$parameter[i],
                comparison$true_value[i],
                comparison$mixture_mean[i],
                comparison$mixture_sd[i],
                comparison$dw_mean[i],
                comparison$dw_sd[i],
                comparison$diff[i]))
}

# ==============================================================================
# Compare diagnostics
# ==============================================================================

cat("\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("DIAGNOSTICS\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")

diag_mixture <- fit_mixture$diagnostic_summary(quiet = TRUE)
diag_dw <- fit_double_woodbury$diagnostic_summary(quiet = TRUE)

cat("\nMixture model:\n")
cat("  Divergences:", sum(diag_mixture$num_divergent), "\n")
cat("  Max treedepth:", sum(diag_mixture$num_max_treedepth), "\n")

cat("\nDouble-Woodbury model:\n")
cat("  Divergences:", sum(diag_dw$num_divergent), "\n")
cat("  Max treedepth:", sum(diag_dw$num_max_treedepth), "\n")

# Rhat comparison
cat("\nMax Rhat (mixture):", max(summary_mixture$rhat, na.rm = TRUE) |> round(3), "\n")
cat("Max Rhat (double-Woodbury):", max(summary_dw$rhat, na.rm = TRUE) |> round(3), "\n")

cat("\nMin ESS bulk (mixture):", min(summary_mixture$ess_bulk, na.rm = TRUE) |> round(0), "\n")
cat("Min ESS bulk (double-Woodbury):", min(summary_dw$ess_bulk, na.rm = TRUE) |> round(0), "\n")
