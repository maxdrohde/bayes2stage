# Comparison test: Woodbury identity vs old Cholesky update
# Compares timing, ESS, and parameter estimates

library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)

# Output file for results - save to project directory
output_file <- file.path(getwd(), "woodbury_comparison_results.txt")
cat("Results will be written to:", output_file, "\n")

write_output <- function(..., append = TRUE) {
    cat(..., "\n", file = output_file, append = append)
    cat(..., "\n")
}

write_output("=" |> rep(70) |> paste(collapse = ""), append = FALSE)
write_output("Woodbury Identity vs Cholesky Update Comparison")
write_output("=" |> rep(70) |> paste(collapse = ""))
write_output(paste("Timestamp:", Sys.time()))
write_output("")

# ==============================================================================
# 1. Generate data
# ==============================================================================
write_output("1. DATA GENERATION")
write_output("-" |> rep(70) |> paste(collapse = ""))

set.seed(777)

N <- 300L
M <- 5L

# True parameters
true_params <- list(
    alpha_main = 1,
    beta_x = 1,
    beta_z = 1,
    beta_t = 2,
    beta_t_x_interaction = 0.3,
    error_sd = 4,
    rand_intercept_sd = 3,
    rand_slope_sd = 1,
    rand_eff_corr = 0.3,
    gamma0 = 1,
    gamma1 = 1,
    gamma_sd = 2
)

write_output(paste("N (subjects):", N))
write_output(paste("M (observations per subject):", M))
write_output(paste("Total observations:", N * M))

# Generate data using bayes2stage
data <- bayes2stage::generate_data(
    N = N,
    M = M,
    alpha_main = true_params$alpha_main,
    beta_x = true_params$beta_x,
    beta_z = true_params$beta_z,
    beta_t = true_params$beta_t,
    beta_t_x_interaction = true_params$beta_t_x_interaction,
    error_sd = true_params$error_sd,
    rand_intercept_sd = true_params$rand_intercept_sd,
    rand_slope_sd = true_params$rand_slope_sd,
    rand_eff_corr = true_params$rand_eff_corr,
    gamma0 = true_params$gamma0,
    gamma1 = true_params$gamma1,
    gamma_sd = true_params$gamma_sd
)

# Apply ODS design
data_ods <- bayes2stage::ods_design(
    data = data,
    sampling_type = "slope",
    cutoff_high = 0.75,
    cutoff_low = 0.25,
    n_sampled = 100,
    prop_high = 0.4,
    prop_middle = 0.2,
    prop_low = 0.4
)

n_obs <- sum(!is.na(data_ods$x[!duplicated(data_ods$id)]))
n_mis <- sum(is.na(data_ods$x[!duplicated(data_ods$id)]))
write_output(paste("Subjects with observed x:", n_obs))
write_output(paste("Subjects with missing x:", n_mis))
write_output("")

# ==============================================================================
# 2. Format data for Stan
# ==============================================================================
data_list <- bayes2stage::format_data_mcmc(
    data = data_ods,
    main_model_formula = "~ z",
    imputation_model_formula = "~ z",
    imputation_distribution = "normal"
)

# ==============================================================================
# 3. Compile models
# ==============================================================================
write_output("2. MODEL COMPILATION")
write_output("-" |> rep(70) |> paste(collapse = ""))

stan_dir <- file.path(getwd(), "src", "stan")

write_output("Compiling NEW model (Woodbury identity)...")
mod_new <- cmdstan_model(file.path(stan_dir, "mixed_effects_imputation_normal_marginalized.stan"))
write_output("  Done.")

write_output("Compiling OLD model (Cholesky update)...")
mod_old <- cmdstan_model(file.path(stan_dir, "mixed_effects_imputation_normal_marginalized_old.stan"))
write_output("  Done.")
write_output("")

# ==============================================================================
# 4. Fit models
# ==============================================================================
write_output("3. MODEL FITTING")
write_output("-" |> rep(70) |> paste(collapse = ""))

n_chains <- 4L
iter_warmup <- 1000L
iter_sampling <- 2000L
seed <- 777L

write_output(paste("Chains:", n_chains))
write_output(paste("Warmup iterations:", iter_warmup))
write_output(paste("Sampling iterations:", iter_sampling))
write_output("")

# Fit OLD model
write_output("Fitting OLD model (Cholesky update)...")
time_old <- system.time({
    fit_old <- mod_old$sample(
        data = data_list,
        seed = seed,
        chains = n_chains,
        parallel_chains = n_chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        refresh = 0,
        show_messages = FALSE
    )
})
write_output(paste("  Elapsed time:", round(time_old["elapsed"], 2), "seconds"))

# Fit NEW model
write_output("Fitting NEW model (Woodbury identity)...")
time_new <- system.time({
    fit_new <- mod_new$sample(
        data = data_list,
        seed = seed,
        chains = n_chains,
        parallel_chains = n_chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        refresh = 0,
        show_messages = FALSE
    )
})
write_output(paste("  Elapsed time:", round(time_new["elapsed"], 2), "seconds"))
write_output("")

# ==============================================================================
# 5. Compare timing
# ==============================================================================
write_output("4. TIMING COMPARISON")
write_output("-" |> rep(70) |> paste(collapse = ""))

speedup <- time_old["elapsed"] / time_new["elapsed"]
write_output(paste("OLD model (Cholesky):", round(time_old["elapsed"], 2), "seconds"))
write_output(paste("NEW model (Woodbury):", round(time_new["elapsed"], 2), "seconds"))
write_output(paste("Speedup factor:", round(speedup, 2), "x"))
write_output("")

# ==============================================================================
# 6. Compare ESS
# ==============================================================================
write_output("5. ESS COMPARISON")
write_output("-" |> rep(70) |> paste(collapse = ""))

params_of_interest <- c("alpha_main", "beta_t", "beta_x", "beta_x_t_interaction",
                        "sigma_main", "alpha_imputation", "sigma_imputation",
                        "sigma_re[1]", "sigma_re[2]")

ess_old <- fit_old$summary(variables = params_of_interest, "ess_bulk", "ess_tail")
ess_new <- fit_new$summary(variables = params_of_interest, "ess_bulk", "ess_tail")

ess_comparison <- ess_old |>
    select(variable, ess_bulk_old = ess_bulk, ess_tail_old = ess_tail) |>
    left_join(
        ess_new |> select(variable, ess_bulk_new = ess_bulk, ess_tail_new = ess_tail),
        by = "variable"
    ) |>
    mutate(
        bulk_ratio = ess_bulk_new / ess_bulk_old,
        tail_ratio = ess_tail_new / ess_tail_old
    )

write_output("ESS Bulk:")
write_output(sprintf("  %-25s %10s %10s %10s", "Parameter", "OLD", "NEW", "Ratio"))
for (i in seq_len(nrow(ess_comparison))) {
    write_output(sprintf("  %-25s %10.0f %10.0f %10.2f",
                         ess_comparison$variable[i],
                         ess_comparison$ess_bulk_old[i],
                         ess_comparison$ess_bulk_new[i],
                         ess_comparison$bulk_ratio[i]))
}

write_output("")
write_output("ESS Tail:")
write_output(sprintf("  %-25s %10s %10s %10s", "Parameter", "OLD", "NEW", "Ratio"))
for (i in seq_len(nrow(ess_comparison))) {
    write_output(sprintf("  %-25s %10.0f %10.0f %10.2f",
                         ess_comparison$variable[i],
                         ess_comparison$ess_tail_old[i],
                         ess_comparison$ess_tail_new[i],
                         ess_comparison$tail_ratio[i]))
}

write_output("")
write_output(paste("Mean ESS bulk ratio (NEW/OLD):", round(mean(ess_comparison$bulk_ratio), 3)))
write_output(paste("Mean ESS tail ratio (NEW/OLD):", round(mean(ess_comparison$tail_ratio), 3)))
write_output("")

# ESS per second
ess_per_sec_old <- mean(ess_comparison$ess_bulk_old) / time_old["elapsed"]
ess_per_sec_new <- mean(ess_comparison$ess_bulk_new) / time_new["elapsed"]
write_output(paste("ESS/second (OLD):", round(ess_per_sec_old, 1)))
write_output(paste("ESS/second (NEW):", round(ess_per_sec_new, 1)))
write_output(paste("ESS/second improvement:", round(ess_per_sec_new / ess_per_sec_old, 2), "x"))
write_output("")

# ==============================================================================
# 7. Compare parameter estimates
# ==============================================================================
write_output("6. PARAMETER ESTIMATES COMPARISON")
write_output("-" |> rep(70) |> paste(collapse = ""))

summary_old <- fit_old$summary(variables = params_of_interest, "mean", "sd", ~quantile(.x, probs = c(0.05, 0.95)))
summary_new <- fit_new$summary(variables = params_of_interest, "mean", "sd", ~quantile(.x, probs = c(0.05, 0.95)))

# Rename columns for easier access
names(summary_old)[names(summary_old) == "5%"] <- "q5"
names(summary_old)[names(summary_old) == "95%"] <- "q95"
names(summary_new)[names(summary_new) == "5%"] <- "q5"
names(summary_new)[names(summary_new) == "95%"] <- "q95"

param_comparison <- summary_old |>
    select(variable, mean_old = mean, sd_old = sd) |>
    left_join(
        summary_new |> select(variable, mean_new = mean, sd_new = sd),
        by = "variable"
    ) |>
    mutate(
        mean_diff = mean_new - mean_old,
        mean_diff_pct = 100 * (mean_new - mean_old) / abs(mean_old),
        sd_ratio = sd_new / sd_old
    )

write_output("Parameter estimates:")
write_output(sprintf("  %-25s %10s %10s %10s %10s", "Parameter", "OLD mean", "NEW mean", "Diff", "Diff %"))
for (i in seq_len(nrow(param_comparison))) {
    write_output(sprintf("  %-25s %10.4f %10.4f %10.4f %9.2f%%",
                         param_comparison$variable[i],
                         param_comparison$mean_old[i],
                         param_comparison$mean_new[i],
                         param_comparison$mean_diff[i],
                         param_comparison$mean_diff_pct[i]))
}

write_output("")
write_output("Posterior SD comparison:")
write_output(sprintf("  %-25s %10s %10s %10s", "Parameter", "OLD SD", "NEW SD", "Ratio"))
for (i in seq_len(nrow(param_comparison))) {
    write_output(sprintf("  %-25s %10.4f %10.4f %10.3f",
                         param_comparison$variable[i],
                         param_comparison$sd_old[i],
                         param_comparison$sd_new[i],
                         param_comparison$sd_ratio[i]))
}

write_output("")
write_output(paste("Max absolute mean difference:", round(max(abs(param_comparison$mean_diff)), 6)))
write_output(paste("Max mean difference (%):", round(max(abs(param_comparison$mean_diff_pct)), 3), "%"))
write_output(paste("Mean SD ratio:", round(mean(param_comparison$sd_ratio), 4)))
write_output("")

# ==============================================================================
# 8. Compare to true values
# ==============================================================================
write_output("7. COMPARISON TO TRUE VALUES")
write_output("-" |> rep(70) |> paste(collapse = ""))

true_values <- data.frame(
    variable = c("alpha_main", "beta_t", "beta_x", "beta_x_t_interaction",
                 "sigma_main", "alpha_imputation", "sigma_imputation",
                 "sigma_re[1]", "sigma_re[2]"),
    true = c(true_params$alpha_main, true_params$beta_t, true_params$beta_x,
             true_params$beta_t_x_interaction, true_params$error_sd,
             true_params$gamma0, true_params$gamma_sd,
             true_params$rand_intercept_sd, true_params$rand_slope_sd)
)

comparison_to_true <- true_values |>
    left_join(summary_old |> select(variable, mean_old = mean, q5_old = q5, q95_old = q95), by = "variable") |>
    left_join(summary_new |> select(variable, mean_new = mean, q5_new = q5, q95_new = q95), by = "variable") |>
    mutate(
        bias_old = mean_old - true,
        bias_new = mean_new - true,
        covers_old = (q5_old <= true) & (true <= q95_old),
        covers_new = (q5_new <= true) & (true <= q95_new)
    )

write_output(sprintf("  %-25s %8s %10s %10s %8s %8s", "Parameter", "True", "OLD mean", "NEW mean", "Cov OLD", "Cov NEW"))
for (i in seq_len(nrow(comparison_to_true))) {
    write_output(sprintf("  %-25s %8.2f %10.4f %10.4f %8s %8s",
                         comparison_to_true$variable[i],
                         comparison_to_true$true[i],
                         comparison_to_true$mean_old[i],
                         comparison_to_true$mean_new[i],
                         ifelse(comparison_to_true$covers_old[i], "Yes", "No"),
                         ifelse(comparison_to_true$covers_new[i], "Yes", "No")))
}

write_output("")
write_output(paste("Coverage (90% CI) - OLD:", sum(comparison_to_true$covers_old), "/", nrow(comparison_to_true)))
write_output(paste("Coverage (90% CI) - NEW:", sum(comparison_to_true$covers_new), "/", nrow(comparison_to_true)))
write_output("")

# ==============================================================================
# 9. Summary
# ==============================================================================
write_output("=" |> rep(70) |> paste(collapse = ""))
write_output("SUMMARY")
write_output("=" |> rep(70) |> paste(collapse = ""))
write_output("")
write_output(paste("Timing speedup:", round(speedup, 2), "x faster with Woodbury"))
write_output(paste("ESS/second improvement:", round(ess_per_sec_new / ess_per_sec_old, 2), "x"))
write_output(paste("Max parameter estimate difference:", round(max(abs(param_comparison$mean_diff_pct)), 4), "%"))
write_output(paste("Models produce equivalent results:", max(abs(param_comparison$mean_diff_pct)) < 1))
write_output("")
write_output(paste("Full results saved to:", output_file))
write_output("")

cat("\n\nDone! Results saved to:", output_file, "\n")
