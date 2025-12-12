#!/usr/bin/env Rscript
################################################################################
# Test Script: Compare Marginal vs Centered Parameterizations
#
# This script tests whether the marginal model improves MCMC mixing for
# the N=2000 case compared to the centered parameterization.
################################################################################

# Suppress package startup messages
suppressPackageStartupMessages({
  library(bayes2stage)
  library(dplyr)
  library(tibble)
})

cat("\n")
cat("=" , rep("=", 78), "\n", sep = "")
cat("Testing MCMC Mixing: Marginal vs Centered Parameterization\n")
cat("=" , rep("=", 78), "\n", sep = "")

################################################################################
# Configuration
################################################################################

# Use smaller N for faster testing, but still large enough to show the issue
TEST_N <- 500L
TEST_M <- 5L
TEST_SAMPLING_FRACTION <- 0.25
TEST_ITER <- 1000L
TEST_WARMUP <- 500L
TEST_CHAINS <- 2L

cat("\nTest Configuration:\n")
cat("  N (subjects):", TEST_N, "\n")
cat("  M (timepoints):", TEST_M, "\n")
cat("  Total observations:", TEST_N * TEST_M, "\n")
cat("  Sampling fraction:", TEST_SAMPLING_FRACTION, "\n")
cat("  Subjects with missing x:", round(TEST_N * (1 - TEST_SAMPLING_FRACTION)), "\n")
cat("  MCMC iterations:", TEST_ITER, "\n")
cat("  MCMC warmup:", TEST_WARMUP, "\n")
cat("  MCMC chains:", TEST_CHAINS, "\n")

################################################################################
# Generate Test Data
################################################################################

cat("\n", rep("-", 80), "\n", sep = "")
cat("Generating test data...\n")

set.seed(12345)

test_data <- bayes2stage::generate_data(
  N = TEST_N,
  M = TEST_M,
  alpha_main = 0,
  beta_x = -1,
  beta_z = -2,
  beta_t = -1,
  beta_t_x_interaction = -0.5,
  beta_t_z_interaction = 0,
  error_sd = 1,
  x_dist = "normal",
  rand_intercept_sd = 4,
  rand_slope_sd = 1,
  rand_eff_corr = 0,
  gamma0 = 0,
  gamma1 = 0.5,
  gamma2 = 0,
  gamma_sd = 1
)

# Apply SRS design to create missing x values
n_sampled <- as.integer(TEST_SAMPLING_FRACTION * TEST_N)
test_data_srs <- bayes2stage::srs_design(test_data, n_sampled)

cat("  Data generated successfully\n")
cat("  Subjects with observed x:", sum(!is.na(test_data_srs$x[!duplicated(test_data_srs$id)])), "\n")
cat("  Subjects with missing x:", sum(is.na(test_data_srs$x[!duplicated(test_data_srs$id)])), "\n")

################################################################################
# Helper function to extract diagnostics
################################################################################

extract_diagnostics <- function(fit, model_name) {
  # Get summary
  summ <- fit$summary()

  # Get main parameters only (exclude random effects and x_mis)
  main_params <- c("alpha_main", "beta_t", "beta_x", "beta_x_t_interaction",
                   "beta[1]", "sigma_main", "sigma_re[1]", "sigma_re[2]",
                   "alpha_imputation", "gamma[1]", "sigma_imputation")

  summ_main <- summ %>%
    filter(variable %in% main_params)

  # Get diagnostics
  diag <- fit$diagnostic_summary()

  tibble(
    model = model_name,
    min_ess_bulk = min(summ_main$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(summ_main$ess_tail, na.rm = TRUE),
    max_rhat = max(summ_main$rhat, na.rm = TRUE),
    n_divergent = sum(diag$num_divergent),
    n_max_treedepth = sum(diag$num_max_treedepth),
    min_ebfmi = min(diag$ebfmi),
    time_seconds = fit$time()$total
  )
}

################################################################################
# Test Centered Parameterization (Original)
################################################################################

cat("\n", rep("-", 80), "\n", sep = "")
cat("Fitting CENTERED model (original parameterization)...\n")
cat("  Expected parameters: ~", 2 * TEST_N + round(TEST_N * (1 - TEST_SAMPLING_FRACTION)) + 12, "\n")

tryCatch({
  t_start <- Sys.time()

  fit_centered <- fit_stan_model(
    data = test_data_srs,
    main_model_covariates = c("z"),
    imputation_model_covariates = c("z"),
    imputation_distribution = "normal",
    parameterization = "centered",
    n_chains = TEST_CHAINS,
    iter_warmup = TEST_WARMUP,
    iter_sampling = TEST_ITER,
    adapt_delta = 0.8,
    parallel_chains = TEST_CHAINS
  )

  t_end <- Sys.time()
  cat("  Completed in", round(difftime(t_end, t_start, units = "secs"), 1), "seconds\n")

  diag_centered <- extract_diagnostics(fit_centered, "centered")

}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  diag_centered <<- tibble(
    model = "centered",
    min_ess_bulk = NA, min_ess_tail = NA, max_rhat = NA,
    n_divergent = NA, n_max_treedepth = NA, min_ebfmi = NA, time_seconds = NA
  )
})

################################################################################
# Test Marginal Parameterization (New)
################################################################################

cat("\n", rep("-", 80), "\n", sep = "")
cat("Fitting MARGINAL model (random effects integrated out)...\n")
cat("  Expected parameters: ~", round(TEST_N * (1 - TEST_SAMPLING_FRACTION)) + 12, "\n")

tryCatch({
  t_start <- Sys.time()

  fit_marginal <- fit_stan_model(
    data = test_data_srs,
    main_model_covariates = c("z"),
    imputation_model_covariates = c("z"),
    imputation_distribution = "normal",
    parameterization = "marginal",
    n_chains = TEST_CHAINS,
    iter_warmup = TEST_WARMUP,
    iter_sampling = TEST_ITER,
    adapt_delta = 0.95,
    max_treedepth = 12,
    parallel_chains = TEST_CHAINS,
    prior_beta_sd = 5.0,
    prior_sigma_rate = 1.0,
    prior_sigma_re_rate = 0.25
  )

  t_end <- Sys.time()
  cat("  Completed in", round(difftime(t_end, t_start, units = "secs"), 1), "seconds\n")

  diag_marginal <- extract_diagnostics(fit_marginal, "marginal")

}, error = function(e) {
  cat("  ERROR:", e$message, "\n")
  diag_marginal <<- tibble(
    model = "marginal",
    min_ess_bulk = NA, min_ess_tail = NA, max_rhat = NA,
    n_divergent = NA, n_max_treedepth = NA, min_ebfmi = NA, time_seconds = NA
  )
})

################################################################################
# Compare Results
################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("RESULTS COMPARISON\n")
cat(rep("=", 80), "\n", sep = "")

results <- bind_rows(diag_centered, diag_marginal)
print(results, width = Inf)

cat("\n")
cat("Interpretation:\n")
cat("  - Higher ESS = better mixing\n")
cat("  - Rhat closer to 1.0 = better convergence (should be < 1.01)\n")
cat("  - Fewer divergences = healthier posterior geometry\n")
cat("  - Higher EBFMI = better energy transitions (should be > 0.3)\n")

cat("\n")
cat(rep("=", 80), "\n", sep = "")

################################################################################
# Parameter Estimates Comparison (if both succeeded)
################################################################################

if (exists("fit_centered") && exists("fit_marginal")) {
  cat("\nParameter Estimates (True values in parentheses):\n")
  cat(rep("-", 80), "\n", sep = "")

  true_vals <- c(
    "alpha_main" = 0,
    "beta_t" = -1,
    "beta_x" = -1,
    "beta_x_t_interaction" = -0.5,
    "beta[1]" = -2,  # beta_z
    "sigma_main" = 1,
    "sigma_re[1]" = 4,
    "sigma_re[2]" = 1
  )

  summ_c <- fit_centered$summary() %>% filter(variable %in% names(true_vals))
  summ_m <- fit_marginal$summary() %>% filter(variable %in% names(true_vals))

  comparison <- summ_c %>%
    select(variable, mean_centered = mean, sd_centered = sd) %>%
    left_join(
      summ_m %>% select(variable, mean_marginal = mean, sd_marginal = sd),
      by = "variable"
    ) %>%
    mutate(true_value = true_vals[variable])

  print(comparison, width = Inf)
}

cat("\nTest completed.\n")
