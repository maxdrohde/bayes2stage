# Validation script for negative-binomial imputation model
# Tests parameter recovery with unbounded count data

library(bayes2stage)
library(dplyr)

# Configuration
N_SUBJECTS <- 300
N_TIMEPOINTS <- 5
SAMPLING_FRACTION <- 0.6
N_CHAINS <- 4
ITER_WARMUP <- 2000
ITER_SAMPLING <- 2000
SEED <- 67890

# True parameter values
TRUE_PARAMS <- list(
  alpha_main = 2.0,
  beta_x = 0.3,  # Effect per unit increase in count
  beta_z = 0.8,
  beta_t = 1.2,
  beta_t_x_interaction = 0.15,
  error_sd = 3.0,
  rand_intercept_sd = 2.5,
  rand_slope_sd = 0.6,
  rand_eff_corr = 0.25,
  gamma0 = 1.0,   # Intercept on log scale (exp(1) = 2.7 mean)
  gamma1 = 0.5,   # Effect of z on log scale
  phi = 3.0       # Dispersion parameter
)

cat("\n")
cat(strrep("=", 80), "\n")
cat("VALIDATING NEGATIVE-BINOMIAL IMPUTATION MODEL\n")
cat(strrep("=", 80), "\n\n")

set.seed(SEED)

# Generate data with negative-binomial count x
df_negbin <- generate_data(
  N = N_SUBJECTS,
  M = N_TIMEPOINTS,
  x_dist = "negative_binomial",
  x_disp_param = TRUE_PARAMS$phi,
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
  gamma1 = TRUE_PARAMS$gamma1
)

cat("Data generated:", N_SUBJECTS, "subjects x", N_TIMEPOINTS, "timepoints\n")
cat("Count range:", min(df_negbin$x), "-",
    max(df_negbin$x[!duplicated(df_negbin$id)]), "\n")
cat("Mean count:", round(mean(df_negbin$x[!duplicated(df_negbin$id)]), 2), "\n\n")

# Apply SRS sampling
sampling_n <- round(N_SUBJECTS * SAMPLING_FRACTION)
df_sampled <- srs_design(df_negbin, n_sampled = sampling_n)

cat("Subjects with observed x:", sampling_n, "\n")
cat("Subjects with missing x:", N_SUBJECTS - sampling_n, "\n\n")

# Fit the model
cat("Fitting negative-binomial imputation model...\n")
cat("Chains:", N_CHAINS, "| Warmup:", ITER_WARMUP, "| Sampling:", ITER_SAMPLING, "\n\n")

fit <- fit_stan_model(
  data = df_sampled,
  main_model_covariates = c("z"),
  imputation_model_covariates = c("z"),
  imputation_distribution = "negative_binomial",
  n_chains = N_CHAINS,
  iter_warmup = ITER_WARMUP,
  iter_sampling = ITER_SAMPLING,
  seed = SEED,
  parallel_chains = N_CHAINS
)

# Check parameter recovery
params_to_check <- list(
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
  phi = TRUE_PARAMS$phi
)

# Summarize each parameter
results <- lapply(names(params_to_check), function(p) {
  draws <- fit$draws(variables = p, format = "matrix")
  true_val <- params_to_check[[p]]
  post_mean <- mean(draws)
  ci <- unname(quantile(draws, c(0.025, 0.975)))
  covered <- true_val >= ci[1] && true_val <= ci[2]

  data.frame(
    parameter = p,
    true = true_val,
    mean = post_mean,
    ci_lower = ci[1],
    ci_upper = ci[2],
    bias = post_mean - true_val,
    covered = covered
  )
}) |> bind_rows()

# Print results
cat("\n")
cat(strrep("=", 80), "\n")
cat("RESULTS\n")
cat(strrep("=", 80), "\n\n")

print(results |>
        mutate(across(where(is.numeric), ~round(.x, 3))) |>
        as.data.frame(),
      row.names = FALSE)

cat("\n")
cat("Coverage rate:", mean(results$covered) * 100, "%\n")
cat("Mean absolute bias:", round(mean(abs(results$bias)), 3), "\n\n")

# Diagnostics
cat("Diagnostics:\n")
print(fit$diagnostic_summary())

# Validation check
if (mean(results$covered) >= 0.80) {
  cat("\nVALIDATION PASSED: Model shows good parameter recovery.\n")
} else {
  cat("\nVALIDATION WARNING: Coverage below 80%. Check model specification.\n")
}
