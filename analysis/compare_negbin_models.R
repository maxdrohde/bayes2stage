# Compare original vs optimized negative binomial Stan models
# Tests that both produce equivalent results and compares timing

library(bayes2stage)
library(cmdstanr)
library(dplyr)

# =============================================================================
# Configuration
# =============================================================================

N_SUBJECTS <- 500
M_TIMEPOINTS <- 5
SAMPLING_FRACTION <- 0.6
N_CHAINS <- 4
ITER_WARMUP <- 5000
ITER_SAMPLING <- 5000
SEED <- 12345

# True parameters
TRUE_PARAMS <- list(
  alpha_main = 2.0,
  beta_x = 0.3,
  beta_z = 0.8,
  beta_t = 1.2,
  beta_x_t_interaction = 0.15,
  error_sd = 3.0,
  rand_intercept_sd = 2.5,
  rand_slope_sd = 0.6,
  gamma0 = 1.0,
  gamma1 = 0.5,
  phi = 3.0
)

# =============================================================================
# Generate Data
# =============================================================================

cat("\n", strrep("=", 70), "\n")
cat("COMPARING NEGATIVE BINOMIAL MODELS\n")
cat(strrep("=", 70), "\n\n")

set.seed(SEED)

df <- generate_data(
  N = N_SUBJECTS,
  M = M_TIMEPOINTS,
  x_dist = "negative_binomial",
  x_disp_param = TRUE_PARAMS$phi,
  alpha_main = TRUE_PARAMS$alpha_main,
  beta_x = TRUE_PARAMS$beta_x,
  beta_z = TRUE_PARAMS$beta_z,
  beta_t = TRUE_PARAMS$beta_t,
  beta_t_x_interaction = TRUE_PARAMS$beta_x_t_interaction,
  error_sd = TRUE_PARAMS$error_sd,
  rand_intercept_sd = TRUE_PARAMS$rand_intercept_sd,
  rand_slope_sd = TRUE_PARAMS$rand_slope_sd,
  rand_eff_corr = 0.25,
  gamma0 = TRUE_PARAMS$gamma0,
  gamma1 = TRUE_PARAMS$gamma1
)

# Apply SRS sampling
n_obs <- round(N_SUBJECTS * SAMPLING_FRACTION)
df_sampled <- srs_design(df, n_sampled = n_obs)

cat("Data generated:\n")
cat("  Subjects:", N_SUBJECTS, "\n")
cat("  Timepoints:", M_TIMEPOINTS, "\n")
cat("  Observed x:", n_obs, "\n")
cat("  Missing x:", N_SUBJECTS - n_obs, "\n")
cat("  Count range:", min(df$x), "-", max(df$x), "\n\n")

# =============================================================================
# Prepare Stan Data
# =============================================================================

# Use format_data_mcmc to prepare data (includes max_x)
stan_data <- format_data_mcmc(
  df_sampled,
  main_model_covariates = "z",
  imputation_model_covariates = "z",
  imputation_distribution = "negative_binomial"
)

cat("Stan data prepared:\n")
cat("  N:", stan_data$N, "\n")
cat("  G:", stan_data$G, "\n")
cat("  G_obs:", stan_data$G_obs, "\n")
cat("  G_mis:", stan_data$G_mis, "\n")
cat("  max_x:", stan_data$max_x, "\n\n")

# =============================================================================
# Compile Both Models
# =============================================================================

cat("Compiling models...\n")

stan_path <- file.path(getwd(), "src/stan")

# Original model
mod_original <- cmdstan_model(
  file.path(stan_path, "mixed_effects_imputation_negative_binomial.stan")
)

# Optimized model
mod_optimized <- cmdstan_model(
  file.path(stan_path, "mixed_effects_imputation_negative_binomial_optimized.stan")
)

cat("  Both models compiled successfully.\n\n")

# =============================================================================
# Fit Original Model
# =============================================================================

cat(strrep("-", 70), "\n")
cat("FITTING ORIGINAL MODEL\n")
cat(strrep("-", 70), "\n\n")

time_original <- system.time({
  fit_original <- mod_original$sample(
    data = stan_data,
    seed = SEED,
    chains = N_CHAINS,
    parallel_chains = N_CHAINS,
    iter_warmup = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    refresh = 0,
    show_messages = FALSE
  )
})

cat("Original model time:", round(time_original["elapsed"], 1), "seconds\n\n")

# =============================================================================
# Fit Optimized Model
# =============================================================================

cat(strrep("-", 70), "\n")
cat("FITTING OPTIMIZED MODEL\n")
cat(strrep("-", 70), "\n\n")

time_optimized <- system.time({
  fit_optimized <- mod_optimized$sample(
    data = stan_data,
    seed = SEED,
    chains = N_CHAINS,
    parallel_chains = N_CHAINS,
    iter_warmup = ITER_WARMUP,
    iter_sampling = ITER_SAMPLING,
    refresh = 0,
    show_messages = FALSE
  )
})

cat("Optimized model time:", round(time_optimized["elapsed"], 1), "seconds\n\n")

# =============================================================================
# Compare Results
# =============================================================================

cat(strrep("=", 70), "\n")
cat("RESULTS COMPARISON\n")
cat(strrep("=", 70), "\n\n")

# Parameters to compare
params <- c("alpha_main", "beta_x", "beta_t", "beta_x_t_interaction",
            "beta[1]", "sigma_main", "sigma_re[1]", "sigma_re[2]",
            "alpha_imputation", "gamma[1]", "phi")

# Extract summaries
get_summary <- function(fit, params) {
  fit$summary(variables = params) |>
    select(variable, mean, sd, q5, q95, rhat, ess_bulk)
}

summary_original <- get_summary(fit_original, params)
summary_optimized <- get_summary(fit_optimized, params)

# Combine for comparison
comparison <- summary_original |>
  rename(mean_orig = mean, sd_orig = sd, rhat_orig = rhat) |>
  select(variable, mean_orig, sd_orig, rhat_orig) |>
  left_join(
    summary_optimized |>
      rename(mean_opt = mean, sd_opt = sd, rhat_opt = rhat) |>
      select(variable, mean_opt, sd_opt, rhat_opt),
    by = "variable"
  ) |>
  mutate(
    mean_diff = mean_opt - mean_orig,
    sd_diff = sd_opt - sd_orig,
    pct_diff = abs(mean_diff / mean_orig) * 100
  )

cat("Parameter Estimates:\n\n")
print(comparison |> mutate(across(where(is.numeric), ~round(.x, 4))))

# =============================================================================
# Timing Summary
# =============================================================================

cat("\n")
cat(strrep("=", 70), "\n")
cat("TIMING SUMMARY\n")
cat(strrep("=", 70), "\n\n")

speedup <- time_original["elapsed"] / time_optimized["elapsed"]

cat("Original model:  ", round(time_original["elapsed"], 1), "seconds\n")
cat("Optimized model: ", round(time_optimized["elapsed"], 1), "seconds\n")
cat("Speedup:         ", round(speedup, 2), "x\n\n")

# =============================================================================
# Validation
# =============================================================================

cat(strrep("=", 70), "\n")
cat("VALIDATION\n")
cat(strrep("=", 70), "\n\n")

# Check if estimates are similar (within MCMC error)
max_pct_diff <- max(comparison$pct_diff, na.rm = TRUE)
max_mean_diff <- max(abs(comparison$mean_diff), na.rm = TRUE)

# R-hat check
rhat_ok_orig <- all(comparison$rhat_orig < 1.05, na.rm = TRUE)
rhat_ok_opt <- all(comparison$rhat_opt < 1.05, na.rm = TRUE)

cat("Maximum absolute difference in means:", round(max_mean_diff, 4), "\n")
cat("Maximum percent difference:", round(max_pct_diff, 2), "%\n")
cat("Original R-hat < 1.05:", rhat_ok_orig, "\n")
cat("Optimized R-hat < 1.05:", rhat_ok_opt, "\n\n")

# Pass/fail
if (max_pct_diff < 5 && rhat_ok_orig && rhat_ok_opt) {
  cat("VALIDATION PASSED: Models produce equivalent results.\n")
} else {
  cat("VALIDATION WARNING: Check results for discrepancies.\n")
}

cat("\nComparison complete.\n")
