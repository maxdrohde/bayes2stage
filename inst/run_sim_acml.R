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
output_file <- glue::glue("results_acml/{i}/draws_{j}.parquet")
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
dir_path <- glue::glue("results_acml/{i}")

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
sampling_type <- params[["sampling_type"]]
cutoff_high <- DEFAULT_CUTOFF_HIGH
cutoff_low <- DEFAULT_CUTOFF_LOW
prop_high <- DEFAULT_PROP_HIGH
prop_middle <- DEFAULT_PROP_MIDDLE
prop_low <- DEFAULT_PROP_LOW

################################################################################
# Run simulation
################################################################################

df <- generate_simulation_data(params)

ods_df <- bayes2stage::ods_design(df,
  sampling_type = sampling_type,
  cutoff_high = cutoff_high,
  cutoff_low = cutoff_low,
  n_sampled = n_sampled,
  prop_high = prop_high,
  prop_middle = prop_middle,
  prop_low = prop_low
)

################################################################################
# Fit models
################################################################################

res <- tryCatch(
  {
    invisible(capture.output({
      acml_ods <- bayes2stage::fit_acml_ods(ods_df,
        cutoff_low = cutoff_low,
        cutoff_high = cutoff_high
      )
    }))
    acml_ods
  },
  error = function(e) {
    return(data.frame(message = e$message))
  }
)

################################################################################
# Add metadata to results file
################################################################################

res$N <- N
res$n_sampled <- n_sampled

res$sim_setting <- i
res$sim_iter <- j

res <- dplyr::bind_cols(
  res,
  params[rep(1, nrow(res)), ]
)

################################################################################
# Save results
################################################################################

arrow::write_parquet(res, output_file)
