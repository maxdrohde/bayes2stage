# Validation script for Stan models
# Tests that both normal and bernoulli imputation models can recover true parameters

library(bayes2stage)
library(dplyr)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

N_SUBJECTS <- 4000
N_TIMEPOINTS <- 10
SAMPLING_FRACTION <- 0.6
N_CHAINS <- 10
ITER_WARMUP <- 2000
ITER_SAMPLING <- 2000
SEED <- 12345

# True parameter values
TRUE_PARAMS <- list(
  alpha_main = 2.0,
  beta_x = 1.5,
  beta_z = 0.8,
  beta_t = 1.2,
  beta_t_x_interaction = 0.3,
  error_sd = 3.0,
  rand_intercept_sd = 2.5,
  rand_slope_sd = 0.6,
  rand_eff_corr = 0.25,
  gamma0 = 0.5,
  gamma1 = 1.0,
  gamma_sd = 1.5
)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

summarize_posterior <- function(fit, param_name, true_value) {
  draws <- fit$draws(variables = param_name, format = "matrix")
  post_mean <- mean(draws)
  post_sd <- sd(draws)
  ci_lower <- unname(quantile(draws, 0.025))
  ci_upper <- unname(quantile(draws, 0.975))
  covered <- true_value >= ci_lower && true_value <= ci_upper
  bias <- post_mean - true_value
  rel_bias <- bias / abs(true_value) * 100

  data.frame(
    parameter = param_name,
    true_value = true_value,
    post_mean = post_mean,
    post_sd = post_sd,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    bias = bias,
    rel_bias_pct = rel_bias,
    covered = covered
  )
}

print_results <- function(results_df, model_name) {
  cat("\n")
  cat(strrep("=", 80), "\n")
  cat("Results for:", model_name, "\n")
  cat(strrep("=", 80), "\n\n")

  print(results_df |>
          mutate(across(where(is.numeric), ~round(.x, 3))) |>
          as.data.frame(),
        row.names = FALSE)

  cat("\n")
  cat("Coverage rate:", mean(results_df$covered) * 100, "%\n")
  cat("Mean absolute bias:", mean(abs(results_df$bias)), "\n")
  cat("Mean absolute relative bias:", mean(abs(results_df$rel_bias_pct)), "%\n")
  cat("\n")
}

# -----------------------------------------------------------------------------
# Validate Normal Imputation Model
# -----------------------------------------------------------------------------

cat("\n")
cat(strrep("#", 80), "\n")
cat("# VALIDATING NORMAL IMPUTATION MODEL\n")
cat(strrep("#", 80), "\n")

set.seed(SEED)

# Generate data with continuous x (normal distribution)
df_normal <- generate_data(
  N = N_SUBJECTS,
  M = N_TIMEPOINTS,
  x_dist = "normal",
  alpha_main = TRUE_PARAMS$alpha_main,
  beta_x = TRUE_PARAMS$beta_x,
  beta_z = TRUE_PARAMS$beta_z,
  beta_t = TRUE_PARAMS$beta_t,
  beta_t_x_interaction = TRUE_PARAMS$beta_t_x_interaction,
  error_sd = TRUE_PARAMS$error_sd,
  rand_intercept_sd = TRUE_PARAMS$rand_intercept_sd,
  rand_slope_sd = TRUE_PARAMS$rand_slope_sd,
  rand_eff_corr = TRUE_PARAMS$rand_eff_corr,
  gamma0 = TRUE_PARAMS$gamma0,
  gamma1 = TRUE_PARAMS$gamma1,
  gamma_sd = TRUE_PARAMS$gamma_sd
)

cat("\nData generated:", N_SUBJECTS, "subjects x", N_TIMEPOINTS, "timepoints\n")
cat("Total observations:", nrow(df_normal), "\n")

# Apply SRS sampling to create missingness
sampling_n <- round(N_SUBJECTS * SAMPLING_FRACTION)
df_normal_sampled <- srs_design(df_normal, sampling_N = sampling_n)

cat("Subjects with observed x:", sampling_n, "\n")
cat("Subjects with missing x:", N_SUBJECTS - sampling_n, "\n")

# Fit the normal imputation model
cat("\nFitting normal imputation model...\n")
cat("Chains:", N_CHAINS, "| Warmup:", ITER_WARMUP, "| Sampling:", ITER_SAMPLING, "\n\n")

fit_normal <- fit_stan_model(
  data = df_normal_sampled,
  main_model_covariates = c("z"),
  imputation_model_covariates = c("z"),
  imputation_distribution = "normal",
  nchains = N_CHAINS,
  iter_warmup = ITER_WARMUP,
  iter_sampling = ITER_SAMPLING,
  seed = SEED,
  parallel_chains = N_CHAINS
)

# Summarize results for normal model
params_to_check_normal <- list(
  alpha_main = TRUE_PARAMS$alpha_main,
  beta_x = TRUE_PARAMS$beta_x,
  beta_t = TRUE_PARAMS$beta_t,
  beta_x_t_interaction = TRUE_PARAMS$beta_t_x_interaction,
  `beta[1]` = TRUE_PARAMS$beta_z,
  sigma_main = TRUE_PARAMS$error_sd,
  `sigma_re[1]` = TRUE_PARAMS$rand_intercept_sd,
  `sigma_re[2]` = TRUE_PARAMS$rand_slope_sd,
  alpha_imputation = TRUE_PARAMS$gamma0,
  `gamma[1]` = TRUE_PARAMS$gamma1,
  sigma_imputation = TRUE_PARAMS$gamma_sd
)

results_normal <- bind_rows(
  lapply(names(params_to_check_normal), function(p) {
    summarize_posterior(fit_normal, p, params_to_check_normal[[p]])
  })
)

print_results(results_normal, "Normal Imputation Model")

# Print diagnostics
cat("Diagnostics:\n")
print(fit_normal$diagnostic_summary())

# -----------------------------------------------------------------------------
# Validate Bernoulli Imputation Model
# -----------------------------------------------------------------------------

cat("\n")
cat(strrep("#", 80), "\n")
cat("# VALIDATING BERNOULLI IMPUTATION MODEL\n")
cat(strrep("#", 80), "\n")

set.seed(SEED + 1)

# True parameters for bernoulli model (logit scale for imputation)
TRUE_PARAMS_BERN <- list(
  alpha_main = 2.0,
  beta_x = 2.0,  # Effect of binary x
  beta_z = 0.8,
  beta_t = 1.2,
  beta_t_x_interaction = 0.4,
  error_sd = 3.0,
  rand_intercept_sd = 2.5,
  rand_slope_sd = 0.6,
  rand_eff_corr = 0.25,
  gamma0 = -0.5,  # Intercept on logit scale (~38% baseline probability)
  gamma1 = 1.2    # Effect of z on logit scale
)

# Generate data with binary x (binomial with size=1)
df_bernoulli <- generate_data(
  N = N_SUBJECTS,
  M = N_TIMEPOINTS,
  x_dist = "binomial",
  x_size = 1,
  alpha_main = TRUE_PARAMS_BERN$alpha_main,
  beta_x = TRUE_PARAMS_BERN$beta_x,
  beta_z = TRUE_PARAMS_BERN$beta_z,
  beta_t = TRUE_PARAMS_BERN$beta_t,
  beta_t_x_interaction = TRUE_PARAMS_BERN$beta_t_x_interaction,
  error_sd = TRUE_PARAMS_BERN$error_sd,
  rand_intercept_sd = TRUE_PARAMS_BERN$rand_intercept_sd,
  rand_slope_sd = TRUE_PARAMS_BERN$rand_slope_sd,
  rand_eff_corr = TRUE_PARAMS_BERN$rand_eff_corr,
  gamma0 = TRUE_PARAMS_BERN$gamma0,
  gamma1 = TRUE_PARAMS_BERN$gamma1
)

cat("\nData generated:", N_SUBJECTS, "subjects x", N_TIMEPOINTS, "timepoints\n")
cat("Total observations:", nrow(df_bernoulli), "\n")
cat("Proportion x=1:", mean(df_bernoulli$x[!duplicated(df_bernoulli$id)]), "\n")

# Apply SRS sampling to create missingness
df_bernoulli_sampled <- srs_design(df_bernoulli, sampling_N = sampling_n)

cat("Subjects with observed x:", sampling_n, "\n")
cat("Subjects with missing x:", N_SUBJECTS - sampling_n, "\n")

# Fit the bernoulli imputation model
cat("\nFitting bernoulli imputation model...\n")
cat("Chains:", N_CHAINS, "| Warmup:", ITER_WARMUP, "| Sampling:", ITER_SAMPLING, "\n\n")

fit_bernoulli <- fit_stan_model(
  data = df_bernoulli_sampled,
  main_model_covariates = c("z"),
  imputation_model_covariates = c("z"),
  imputation_distribution = "bernoulli",
  nchains = N_CHAINS,
  iter_warmup = ITER_WARMUP,
  iter_sampling = ITER_SAMPLING,
  seed = SEED,
  parallel_chains = N_CHAINS
)

# Summarize results for bernoulli model
params_to_check_bernoulli <- list(
  alpha_main = TRUE_PARAMS_BERN$alpha_main,
  beta_x = TRUE_PARAMS_BERN$beta_x,
  beta_t = TRUE_PARAMS_BERN$beta_t,
  beta_x_t_interaction = TRUE_PARAMS_BERN$beta_t_x_interaction,
  `beta[1]` = TRUE_PARAMS_BERN$beta_z,
  sigma_main = TRUE_PARAMS_BERN$error_sd,
  `sigma_re[1]` = TRUE_PARAMS_BERN$rand_intercept_sd,
  `sigma_re[2]` = TRUE_PARAMS_BERN$rand_slope_sd,
  alpha_imputation = TRUE_PARAMS_BERN$gamma0,
  `gamma[1]` = TRUE_PARAMS_BERN$gamma1
)

results_bernoulli <- bind_rows(
  lapply(names(params_to_check_bernoulli), function(p) {
    summarize_posterior(fit_bernoulli, p, params_to_check_bernoulli[[p]])
  })
)

print_results(results_bernoulli, "Bernoulli Imputation Model")

# Print diagnostics
cat("Diagnostics:\n")
print(fit_bernoulli$diagnostic_summary())

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

cat("\n")
cat(strrep("#", 80), "\n")
cat("# VALIDATION SUMMARY\n")
cat(strrep("#", 80), "\n\n")

cat("Normal Model:\n")
cat("  - Coverage:", mean(results_normal$covered) * 100, "%\n")
cat("  - Mean |bias|:", round(mean(abs(results_normal$bias)), 3), "\n")
cat("  - Mean |rel bias|:", round(mean(abs(results_normal$rel_bias_pct)), 1), "%\n\n")

cat("Bernoulli Model:\n")
cat("  - Coverage:", mean(results_bernoulli$covered) * 100, "%\n")
cat("  - Mean |bias|:", round(mean(abs(results_bernoulli$bias)), 3), "\n")
cat("  - Mean |rel bias|:", round(mean(abs(results_bernoulli$rel_bias_pct)), 1), "%\n\n")

# Check if validation passed (coverage should be close to 95%)
normal_pass <- mean(results_normal$covered) >= 0.80
bernoulli_pass <- mean(results_bernoulli$covered) >= 0.80

if (normal_pass && bernoulli_pass) {
  cat("VALIDATION PASSED: Both models show acceptable parameter recovery.\n")
} else {
  cat("VALIDATION WARNING: One or more models may have issues with parameter recovery.\n")
  if (!normal_pass) cat("  - Normal model coverage below 80%\n")
  if (!bernoulli_pass) cat("  - Bernoulli model coverage below 80%\n")
}
