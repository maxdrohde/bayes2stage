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
if (fs::file_exists(output_file)) {
  cli::cli_alert_info("Results already exist for setting {i}, iteration {j}")
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

# Using explicit namespacing instead of library() calls

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
if (!fs::dir_exists(dir_path)) {
  fs::dir_create(dir_path)
}

################################################################################
# Simulation scenario grid
################################################################################

grid <- get_simulation_grid()
params <- dplyr::slice(grid, i)

cli::cli_h1("Simulation: i={i}, j={j}, seed={seed}")
cli::cli_h2("Parameters")
cli::cli_dl(as.list(params))
cli::cli_rule()


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
  fit_args <- c(
    list(
      data = data,
      main_model_covariates = DEFAULT_MAIN_COVARIATES,
      imputation_model_covariates = DEFAULT_IMPUTATION_COVARIATES,
      imputation_distribution = DEFAULT_STAN_DISTRIBUTION,
      parameterization = DEFAULT_PARAMETERIZATION,
      seed = seed
    ),
    DEFAULT_INFERENCE_ARGS
  )
  fit <- do.call(bayes2stage::fit_stan_model, fit_args)
  return(fit)
}

fit_acml <- function(data, cutoff_low, cutoff_high) {
  result <- tryCatch(
    {
      invisible(utils::capture.output({
        acml_fit <- bayes2stage::fit_acml_ods(data,
          cutoff_low = cutoff_low,
          cutoff_high = cutoff_high
        )
      }))
      # Use model_summary() to standardize output format
      out <- bayes2stage::model_summary(acml_fit)
      out$type <- "acml_ods"
      out
    },
    error = \(e) {
      data.frame(
        parameter = NA_character_,
        mean = NA_real_,
        sd = NA_real_,
        type = "acml_ods",
        divergent_transitions = NA_integer_,
        max_treedepth_exceeded = NA_integer_,
        ebfmi_min = NA_real_,
        message = e$message
      )
    }
  )
  return(result)
}

df <- generate_simulation_data(params)

datasets <- list()

if ("full" %in% SAMPLING_DESIGNS) {
  datasets[["full"]] <- df
}

if (any(c("srs", "srs_no_imp") %in% SAMPLING_DESIGNS)) {
  srs_df <- bayes2stage::srs_design(df, n_sampled)
  if ("srs" %in% SAMPLING_DESIGNS) {
    datasets[["srs"]] <- srs_df
  }
  if ("srs_no_imp" %in% SAMPLING_DESIGNS) {
    srs_df_no_imp <- tidyr::drop_na(srs_df)
    srs_df_no_imp$id <- match(
      srs_df_no_imp$id,
      unique(srs_df_no_imp$id)
    )
    datasets[["srs_no_imp"]] <- srs_df_no_imp
  }
}

if ("ods" %in% SAMPLING_DESIGNS || FIT_ACML) {
  ods_df <- bayes2stage::ods_design(df,
    sampling_type = params[["sampling_type"]],
    cutoff_high = DEFAULT_CUTOFF_HIGH,
    cutoff_low = DEFAULT_CUTOFF_LOW,
    n_sampled = n_sampled,
    prop_high = DEFAULT_PROP_HIGH,
    prop_middle = DEFAULT_PROP_MIDDLE,
    prop_low = DEFAULT_PROP_LOW
  )
  if ("ods" %in% SAMPLING_DESIGNS) {
    datasets[["ods"]] <- ods_df
  }
}

if ("bds" %in% SAMPLING_DESIGNS) {
  bds_df <- bayes2stage::bds_design(df,
    fixed_effects_formula = DEFAULT_BDS_FORMULA,
    sampling_type = params[["sampling_type"]],
    cutoff_high = DEFAULT_CUTOFF_HIGH,
    cutoff_low = DEFAULT_CUTOFF_LOW,
    n_sampled = n_sampled,
    prop_high = DEFAULT_PROP_HIGH,
    prop_middle = DEFAULT_PROP_MIDDLE,
    prop_low = DEFAULT_PROP_LOW
  )
  datasets[["bds"]] <- bds_df
}

# Fit each of the Stan models
results_list <- list()
for (type in names(datasets)) {
  m <- fit_model(datasets[[type]])
  result <- bayes2stage::model_summary(m)
  result$type <- type
  results_list[[type]] <- result
}

# Fit ACML ODS (uses different fitting method, not Stan)
if (FIT_ACML) {
  acml_result <- fit_acml(ods_df, DEFAULT_CUTOFF_LOW, DEFAULT_CUTOFF_HIGH)
  results_list[["acml_ods"]] <- acml_result
}

res <- dplyr::bind_rows(results_list)

################################################################################
# Add metadata to results file
################################################################################

res$n_sampled <- n_sampled

res$sim_setting <- i
res$sim_iter <- j

res$inference_method <- DEFAULT_INFERENCE_ARGS$inference_method
res$stan_parameterization <- DEFAULT_PARAMETERIZATION

# Add method-specific config metadata
if (DEFAULT_INFERENCE_ARGS$inference_method == "mcmc") {
  res$MCMC_iterations <- DEFAULT_INFERENCE_ARGS$iter_sampling
  res$MCMC_burnin <- DEFAULT_INFERENCE_ARGS$iter_warmup
  res$MCMC_chains <- DEFAULT_INFERENCE_ARGS$n_chains
  res$MCMC_parallel_chains <- DEFAULT_INFERENCE_ARGS$parallel_chains
  res$MCMC_adapt_delta <- DEFAULT_INFERENCE_ARGS$adapt_delta
  res$stan_use_pathfinder <- DEFAULT_INFERENCE_ARGS$use_pathfinder_init
}

res <- dplyr::bind_cols(
  res,
  params[rep(1, nrow(res)), ]
)

################################################################################
# Save results
################################################################################

arrow::write_parquet(res, output_file)
