library(bayes2stage)
library(tidyverse)
library(glue)
library(cmdstanr)

# Ensure cmdstan is installed and path is set
# cmdstanr::check_cmdstan_toolchain()

source("inst/simulation_config.R")

# --- Setup Data ---
cat("Generating data...\n")
grid <- get_simulation_grid()
# Use the first setting from the grid
params <- dplyr::slice(grid, 1)

set.seed(123)
N <- 2000
params[["N"]] <- N
df <- generate_simulation_data(params)

# Create ODS design (realistic scenario)
n_sampled <- as.integer(0.25 * N)
ods_df <- bayes2stage::ods_design(df,
    sampling_type = "intercept",
    cutoff_high = DEFAULT_CUTOFF_HIGH,
    cutoff_low = DEFAULT_CUTOFF_LOW,
    n_sampled = n_sampled,
    prop_high = DEFAULT_PROP_HIGH,
    prop_middle = DEFAULT_PROP_MIDDLE,
    prop_low = DEFAULT_PROP_LOW
)

cat("Data generated. N =", N, ", n_sampled =", n_sampled, "\n")

# --- Define Benchmark Function ---
run_benchmark <- function(use_pathfinder) {
    cat("\n--------------------------------------------------------\n")
    cat("Running with use_pathfinder_init =", use_pathfinder, "\n")
    cat("--------------------------------------------------------\n")

    start_time <- Sys.time()

    fit <- fit_stan_model(
        data = ods_df,
        main_model_covariates = c("z"),
        imputation_model_covariates = c("z"),
        imputation_distribution = "normal",
        parameterization = "noncentered",
        use_pathfinder_init = use_pathfinder,
        pathfinder_num_paths = 20L,
        n_chains = 4L,
        iter_warmup = 2000L,
        iter_sampling = 2000L,
        adapt_delta = 0.8,
        seed = 12345,
        parallel_chains = 4L
    )

    end_time <- Sys.time()
    total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Extract diagnostics
    summary_stats <- fit$summary()

    min_ess <- min(summary_stats$ess_bulk, na.rm = TRUE)
    max_rhat <- max(summary_stats$rhat, na.rm = TRUE)
    divergences <- sum(fit$sampler_diagnostics()[, , "divergent__"])

    # Get Stan's reported timing
    time_info <- fit$time()
    warmup_time <- sum(time_info$chains$warmup)
    sampling_time <- sum(time_info$chains$sampling)

    cat("\nResults for use_pathfinder_init =", use_pathfinder, ":\n")
    cat("Total Wall Time:  ", round(total_time, 2), "s\n")
    cat("Total Warmup (sum):    ", round(warmup_time, 2), "s\n")
    cat("Total Sampling (sum):  ", round(sampling_time, 2), "s\n")
    cat("Min ESS (bulk):   ", round(min_ess, 0), "\n")
    cat("Max R-hat:        ", round(max_rhat, 4), "\n")
    cat("Divergences:      ", divergences, "\n")

    return(list(
        use_pathfinder = use_pathfinder,
        total_time = total_time,
        warmup_time = warmup_time,
        sampling_time = sampling_time,
        min_ess = min_ess,
        max_rhat = max_rhat,
        divergences = divergences
    ))
}

# --- Run Benchmarks ---
# res_standard <- run_benchmark(use_pathfinder = FALSE)
res_pathfinder <- run_benchmark(use_pathfinder = TRUE)

# --- Summary Comparison ---
cat("\n\n========================================================\n")
cat("BENCHMARK SUMMARY\n")
cat("========================================================\n")
cat("Standard run skipped.\n")
cat("Pathfinder run (20 paths) completed.\n")

print(res_pathfinder)
