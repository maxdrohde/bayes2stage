################################################################################
# Simulation Configuration
#
# Shared configuration for simulation scripts. Source this file to access:
#   - get_simulation_grid(): Returns the parameter grid
#   - generate_simulation_data(params, N): Generates data from parameters
#   - get_true_parameter_values(): Maps grid params to output names
#   - format_parameter_name(): Maps output names to display names
#   - clean_type_names(): Maps raw type names to display names
#   - TYPE_LEVELS: Factor levels for analysis types
#   - Default constants (N, sampling parameters, MCMC settings, etc.)
################################################################################

################################################################################
# Default Constants
################################################################################

DEFAULT_CUTOFF_HIGH <- 0.9
DEFAULT_CUTOFF_LOW <- 0.1
DEFAULT_PROP_HIGH <- 0.40
DEFAULT_PROP_MIDDLE <- 0.20
DEFAULT_PROP_LOW <- 0.40

DEFAULT_MCMC_ITERATIONS <- 5000L
DEFAULT_MCMC_BURNIN <- 5000L
DEFAULT_MCMC_CHAINS <- 4L

DEFAULT_N_BOOT_REPS <- 5000L

SHOW_MCMC_DIAGNOSTICS <- TRUE

DEFAULT_STAN_DISTRIBUTION <- "normal"

################################################################################
# Simulation Grid
################################################################################
#
#' Get the simulation parameter grid
#'
#' @return A tibble with all combinations of simulation parameters
get_simulation_grid <- function() {
  tidyr::crossing(
    sampling_type = c("intercept"),
    beta_x = c(-1),
    beta_z = c(-2),
    beta_t_x_interaction = c(-0.5),
    beta_t_z_interaction = c(0),
    M = c(5L),
    alpha_main = 0,
    beta_t = -1,
    error_sd = 1,
    x_dist = "normal",
    x_size = NULL,
    x_disp_param = NULL,
    rand_intercept_sd = 4,
    rand_slope_sd = 1,
    rand_eff_corr = 0,
    gamma0 = c(0),
    gamma1 = c(0.5),
    gamma2 = c(0),
    gamma_sd = c(1),
    # New parameters added to grid
    N = c(2000L),
    sampling_fraction = c(0.25)
  )
}

################################################################################
# Analysis Constants
################################################################################

TYPE_LEVELS <- c(
  "BDS",
  "ODS",
  "ACML ODS",
  "SRS",
  "SRS\n(no imp)"
)

#' Clean type names for display
#'
#' Maps raw type names to formatted display names
#'
#' @param x Character vector of type names
#' @return Character vector of cleaned type names
clean_type_names <- function(x) {
  dplyr::case_when(
    x == "bds" ~ "BDS",
    x == "ods" ~ "ODS",
    x == "ACML ODS" ~ "ACML ODS",
    x == "srs" ~ "SRS",
    x == "srs_no_imp" ~ "SRS\n(no imp)",
    TRUE ~ x
  )
}

################################################################################
# Parameter Name Mappings
################################################################################

#' Get true parameter values for analysis scripts
#'
#' Maps grid parameters to their output names (e.g., beta_z -> "beta[1]")
#'
#' @param grid_row A single row from the simulation grid (default: row 1)
#' @return A named vector of true parameter values
get_true_parameter_values <- function(grid_row = NULL) {
  if (is.null(grid_row)) {
    grid_row <- dplyr::slice(get_simulation_grid(), 1)
  }

  c(
    "beta_x" = grid_row[["beta_x"]],
    "beta[1]" = grid_row[["beta_z"]],
    "beta_x_t_interaction" = grid_row[["beta_t_x_interaction"]],
    "beta_t" = grid_row[["beta_t"]],
    "alpha_main" = grid_row[["alpha_main"]],
    "sigma_main" = grid_row[["error_sd"]],
    "sigma_re[1]" = grid_row[["rand_intercept_sd"]],
    "sigma_re[2]" = grid_row[["rand_slope_sd"]],
    "alpha_imputation" = grid_row[["gamma0"]],
    "gamma[1]" = grid_row[["gamma1"]]
  )
}

#' Get true parameter values for all simulation settings
#'
#' Creates a tibble with true parameter values for each simulation setting,
#' enabling correct analysis when running multiple settings with different
#' true values.
#'
#' @return A tibble with columns: sim_setting, parameter, true_value
get_all_true_parameter_values <- function() {
  grid <- get_simulation_grid()
  purrr::map_df(seq_len(nrow(grid)), function(i) {
    row <- dplyr::slice(grid, i)
    params <- get_true_parameter_values(row)
    tibble::tibble(
      sim_setting = i,
      parameter = names(params),
      true_value = unname(params)
    )
  })
}

#' Format parameter names for display
#'
#' Maps output parameter names to human-readable display names
#'
#' @param x Character vector of parameter names
#' @return Character vector of formatted display names
format_parameter_name <- function(x) {
  dplyr::case_when(
    x == "beta_x" ~ "Beta x",
    x == "beta_t" ~ "Beta t",
    x == "beta_x_t_interaction" ~ "Beta x:t Interaction",
    x == "alpha_main" ~ "Intercept (Main)",
    x == "sigma_main" ~ "Error SD (Main)",
    x == "beta[1]" ~ "Beta z",
    x == "sigma_re[1]" ~ "Random Intercept SD",
    x == "sigma_re[2]" ~ "Random Slope SD",
    x == "alpha_imputation" ~ "Intercept (Imputation)",
    x == "gamma[1]" ~ "Gamma_1",
    x == "sigma_imputation" ~ "Error SD (Imputation)",
    TRUE ~ x
  )
}

################################################################################
# Data Generation
################################################################################

#' Generate simulation data from parameters
#'
#' @param params A single row from the simulation grid (must include N)
#' @return A data frame of simulated data
generate_simulation_data <- function(params) {
  bayes2stage::generate_data(
    N = params[["N"]],
    M = params[["M"]],
    alpha_main = params[["alpha_main"]],
    beta_x = params[["beta_x"]],
    beta_z = params[["beta_z"]],
    beta_t = params[["beta_t"]],
    beta_t_x_interaction = params[["beta_t_x_interaction"]],
    beta_t_z_interaction = params[["beta_t_z_interaction"]],
    error_sd = params[["error_sd"]],
    x_dist = params[["x_dist"]],
    x_size = params[["x_size"]],
    x_disp_param = params[["x_disp_param"]],
    rand_intercept_sd = params[["rand_intercept_sd"]],
    rand_slope_sd = params[["rand_slope_sd"]],
    rand_eff_corr = params[["rand_eff_corr"]],
    gamma0 = params[["gamma0"]],
    gamma1 = params[["gamma1"]],
    gamma2 = params[["gamma2"]],
    gamma_sd = params[["gamma_sd"]]
  )
}
