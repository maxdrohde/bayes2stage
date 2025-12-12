################################################################################
# Command Line Arguments
################################################################################

# Read in command line arguments
args <- commandArgs(trailingOnly = TRUE)

# i = Simulation scenario
i <- as.integer(args[[1]])

# j = Iteration
j <- as.integer(args[[2]])

# If the results file already exists, quit
output_file <- glue::glue("results/{i}/draws_{j}.parquet")
if (file.exists(output_file)) {
  cat("Results already exist for setting", i, "iteration", j, "\n")
  quit(save = "no")
}

################################################################################
# Seed using Cantor Mapping
################################################################################

seed <- bayes2stage::cantor_seed(i, j)
set.seed(seed)

################################################################################
# Packages
################################################################################

library(bayes2stage)
library(tidyverse)
library(fs)
library(glue)

################################################################################
# Shared Configuration
################################################################################

source("simulation_config.R")

################################################################################
# Setup
################################################################################

# Each simulation scenario gets a unique results folder
dir_path <- glue::glue("results/{i}")

# Create the results_{i} folder if needed
if (!dir_exists(dir_path)) {
  dir_create(dir_path)
}

################################################################################
# Simulation scenario grid
################################################################################

grid <- get_simulation_grid()
params <- dplyr::slice(grid, i)

print("##########################################################################################")
print(glue::glue("i = {i}"))
print(glue::glue("j = {j}"))
print(glue::glue("seed = {seed}"))
print("##########################################################################################")
print("PARAMETERS:")
print("##########################################################################################")
print(params, n = Inf, width = Inf)
print("##########################################################################################")


################################################################################
# Simulation parameters
################################################################################

N <- params[["N"]]
n_sampled <- as.integer(params[["sampling_fraction"]] * N)

# Define the sampling scheme for ODS and BDS
################################################################################
# Run simulation
################################################################################

fit_model <- function(data) {
  fit_stan_model(
    data = data,
    main_model_covariates = c("z"),
    imputation_model_covariates = c("z"),
    imputation_distribution = DEFAULT_STAN_DISTRIBUTION,
    n_chains = DEFAULT_MCMC_CHAINS,
    iter_warmup = DEFAULT_MCMC_BURNIN,
    iter_sampling = DEFAULT_MCMC_ITERATIONS
  )
}

df <- generate_simulation_data(params)

srs_df <- bayes2stage::srs_design(df, n_sampled)

srs_df_no_imp <- tidyr::drop_na(srs_df)

# Renumber IDs to be sequential integers
srs_df_no_imp$id <- match(
  srs_df_no_imp$id,
  unique(srs_df_no_imp$id)
)

ods_df <- bayes2stage::ods_design(df,
  sampling_type = params[["sampling_type"]],
  cutoff_high = DEFAULT_CUTOFF_HIGH,
  cutoff_low = DEFAULT_CUTOFF_LOW,
  n_sampled = n_sampled,
  prop_high = DEFAULT_PROP_HIGH,
  prop_middle = DEFAULT_PROP_MIDDLE,
  prop_low = DEFAULT_PROP_LOW
)

bds_df <- bayes2stage::bds_design(df,
  fixed_effects_formula = y ~ t + z,
  sampling_type = params[["sampling_type"]],
  cutoff_high = DEFAULT_CUTOFF_HIGH,
  cutoff_low = DEFAULT_CUTOFF_LOW,
  n_sampled = n_sampled,
  prop_high = DEFAULT_PROP_HIGH,
  prop_middle = DEFAULT_PROP_MIDDLE,
  prop_low = DEFAULT_PROP_LOW
)

################################################################################
# Fit models
################################################################################

datasets <- list(
  # full = df,
  srs = srs_df,
  srs_no_imp = srs_df_no_imp,
  ods = ods_df
  # bds = bds_df
)

# Fit each of the models
results_list <- list()
for (type in names(datasets)) {
  m <- fit_model(datasets[[type]])
  result <- bayes2stage::mcmc_summary(m, j)
  result$type <- type
  results_list[[type]] <- result
}

res <- dplyr::bind_rows(results_list)

################################################################################
# Add metadata to results file
################################################################################

res$N <- N
res$n_sampled <- n_sampled

res$sim_setting <- i
res$sim_iter <- j

res$MCMC_iterations <- DEFAULT_MCMC_ITERATIONS
res$MCMC_burnin <- DEFAULT_MCMC_BURNIN
res$MCMC_chains <- DEFAULT_MCMC_CHAINS

res <- dplyr::bind_cols(
  res,
  params[rep(1, nrow(res)), ]
)

################################################################################
# Save results
################################################################################

arrow::write_parquet(res, output_file)
