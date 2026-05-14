################################################################################
# Simulation Configuration
#
# Shared configuration for simulation scripts. Source this file to access:
#   - get_simulation_grid(): Returns the parameter grid
#   - generate_simulation_data(params): Generates data from parameters
#   - get_true_parameter_values(): Maps grid params to output names
#   - format_parameter_name(): Maps output names to display names
#   - clean_type_names(): Maps raw type names to display names
#   - TYPE_LEVELS: Factor levels for analysis types
################################################################################

################################################################################
#
#                    SETTINGS YOU'LL COMMONLY CHANGE
#
################################################################################

# ==============================================================================
# Stan Inference Configuration
# ==============================================================================
# Options: "mcmc", "pathfinder", "laplace", "optimize"
#   - mcmc:       Full MCMC (gold standard, slowest)
#   - pathfinder: Fast variational approximation
#   - laplace:    Laplace approximation (very fast, good accuracy)
#   - optimize:   MAP point estimate only (fastest, no uncertainty)

# Arguments passed directly to fit_stan_model()
DEFAULT_INFERENCE_ARGS <- list(
    inference_method = "mcmc",
    n_chains = 4L,
    iter_warmup = 500L,
    iter_sampling = 500L,
    parallel_chains = 1L,
    adapt_delta = 0.8,
    use_pathfinder_init = FALSE
)

# ==============================================================================
# Which Analyses to Run
# ==============================================================================
# Sampling design options: "full", "srs", "srs_no_imp", "ods", "bds"

SAMPLING_DESIGNS <- c("ods", "bds", "srs", "srs_no_imp")

# Whether to fit the ACML model (requires ODS data)
FIT_ACML <- TRUE


################################################################################
#
#                    SETTINGS YOU'LL RARELY CHANGE
#
################################################################################

# ==============================================================================
# ODS/BDS Sampling Cutoffs
# ==============================================================================

DEFAULT_CUTOFF_HIGH <- 0.9
DEFAULT_CUTOFF_LOW <- 0.1
DEFAULT_PROP_HIGH <- 0.40
DEFAULT_PROP_MIDDLE <- 0.20
DEFAULT_PROP_LOW <- 0.40

# ==============================================================================
# Stan Model Specification
# ==============================================================================
# Distribution for x: "normal", "bernoulli", "beta_binomial", "negative_binomial"
# Set DEFAULT_MIXTURE_COMPONENTS to an integer >= 2 to use mixture imputation
# (only for normal distribution), or NULL for single-component (default)

DEFAULT_STAN_DISTRIBUTION <- "bernoulli"
DEFAULT_MIXTURE_COMPONENTS <- NULL

# Covariates for main outcome model and imputation model
DEFAULT_MAIN_FORMULA <- "~ z"
DEFAULT_IMPUTATION_FORMULA <- "~ z"
DEFAULT_BDS_FORMULA <- y ~ t + z

# ==============================================================================
# Analysis/Diagnostics Settings
# ==============================================================================

DEFAULT_CI_LEVEL <- 0.95
DEFAULT_RHAT_THRESHOLD <- 1.01
DEFAULT_ESS_THRESHOLD <- 400L
DEFAULT_N_BOOT_REPS <- 5000L

# Relative efficiency baseline
DEFAULT_BASELINE_TYPE <- "SRS\n(no imp)"

# ==============================================================================
# Plotting Settings
# ==============================================================================

DEFAULT_N_DATA_SAMPLES <- 5L
DEFAULT_PLOT_SUBSET_SIZE <- 200L


################################################################################
#
#                         SIMULATION GRID
#
#   Define the parameter combinations to simulate
#   Each row = one simulation setting
#
################################################################################

#' Get the simulation parameter grid
#'
#' Defaults match Schildcrout 2019 Web Appendix simulation setup.
#'
#' @return A tibble with all combinations of simulation parameters
get_simulation_grid <- function() {
    grid <- tidyr::crossing(
        # Sample size and design
        N = c(300L),
        sampling_fraction = c(0.25),
        sampling_type = c("intercept"),

        # Main model coefficients
        # Paper: (β_0, β_s, β_t, β_st, β_c) = (75, -0.5, -1, -0.5, -2)
        alpha_main = c(75),
        beta_x = c(-0.5),
        # beta_z: 0 (none), -2 (moderate), -4 (strong z effect on intercept)
        beta_z = c(0, -2, -4),
        beta_t = c(-1),
        beta_t_x_interaction = c(-0.5),
        beta_t_z_interaction = c(0),

        # Residual error (fixed at low noise)
        error_sd = c(sqrt(12.25)),

        # Random effects
        # rand_intercept_sd: controls whether ODS selects on z vs random effects
        # Larger values mean ODS selects more on random effects, less on z
        rand_intercept_sd = c(1, 3, 6, 9),
        rand_slope_sd = c(sqrt(1.56)),
        rand_eff_corr = c(0),

        # X distribution: snp ~ Bernoulli(0.25)
        x_dist = c("bernoulli"),
        direction = c("z_given_x"),
        x_prevalence = c(0.25),
        x_size = c(NA_integer_),
        x_disp_param = c(NA_real_),

        # Imputation model: c | snp ~ N(0.25 + gamma1*snp, 1)
        # gamma1 = 0: z independent of x (ACML should not help)
        # gamma1 > 0: z predicts x (ACML should help, more as gamma1 increases)
        gamma0 = c(0.25),
        gamma1 = c(0, 0.5, 1.0, 2.0),
        gamma2 = c(0),
        gamma_sd = c(1),

        # Time scale: integer visit numbers (0, 1, 2, ..., M-1)
        time_scale = c("integer"),

        # M as list column to support variable observations per subject
        # Paper: M sampled from {4, 5, 6}
        M = list(c(4L, 5L, 6L))
    )
    return(grid)
}


################################################################################
#
#                    HELPER FUNCTIONS (Don't Edit)
#
################################################################################

# ==============================================================================
# Type/Parameter Name Formatting
# ==============================================================================

TYPE_LEVELS <- c(
    "Full",
    "BDS",
    "ODS",
    "ACML ODS",
    "SRS",
    "SRS\n(no imp)"
)

#' Clean type names for display
#'
#' @param x Character vector of type names
#' @return Character vector of cleaned type names
clean_type_names <- function(x) {
    result <- dplyr::case_when(
        x == "full" ~ "Full",
        x == "bds" ~ "BDS",
        x == "ods" ~ "ODS",
        x == "acml_ods" ~ "ACML ODS",
        x == "srs" ~ "SRS",
        x == "srs_no_imp" ~ "SRS\n(no imp)",
        TRUE ~ x
    )
    return(result)
}

#' Format parameter names for display
#'
#' @param x Character vector of Stan parameter names
#' @return Character vector of formatted display names
format_parameter_name <- function(x) {
    result <- bayes2stage::translate_parameter_names(
        x,
        from = "stan",
        to = "display"
    )
    return(result)
}

# ==============================================================================
# True Parameter Values (for analysis)
# ==============================================================================

#' Get true parameter values for analysis scripts
#'
#' Maps grid parameters to their output names (e.g., beta_z -> "beta[1]")
#' using the centralized parameter registry.
#'
#' @param grid_row A single row from the simulation grid (default: row 1)
#' @return A named vector of true parameter values (names are Stan parameter names)
get_true_parameter_values <- function(grid_row = NULL) {
    if (is.null(grid_row)) {
        grid_row <- dplyr::slice(get_simulation_grid(), 1)
    }

    registry <- bayes2stage::get_parameter_registry()

    # Extract values from grid_row using grid_name, name with stan_name
    values <- purrr::map_dbl(seq_len(nrow(registry)), \(i) {
        gn <- registry$grid_name[i]
        if (gn %in% names(grid_row)) {
            as.numeric(grid_row[[gn]])
        } else {
            NA_real_
        }
    })
    names(values) <- registry$stan_name
    values <- values[!is.na(values)]
    return(values)
}

#' Get true parameter values for all simulation settings
#'
#' @return A tibble with columns: sim_setting, parameter, true_value
get_all_true_parameter_values <- function() {
    grid <- get_simulation_grid()
    result <- purrr::map_df(seq_len(nrow(grid)), \(i) {
        row <- dplyr::slice(grid, i)
        params <- get_true_parameter_values(row)
        tibble::tibble(
            sim_setting = i,
            parameter = names(params),
            true_value = unname(params)
        )
    })
    return(result)
}

# ==============================================================================
# Data Generation
# ==============================================================================

#' Generate simulation data from parameters
#'
#' @param params A single row from the simulation grid
#' @return A data frame of simulated data
generate_simulation_data <- function(params) {
    data <- bayes2stage::generate_data(
        N = params[["N"]],
        M = params[["M"]][[1]],
        alpha_main = params[["alpha_main"]],
        beta_x = params[["beta_x"]],
        beta_z = params[["beta_z"]],
        beta_t = params[["beta_t"]],
        beta_t_x_interaction = params[["beta_t_x_interaction"]],
        beta_t_z_interaction = params[["beta_t_z_interaction"]],
        error_sd = params[["error_sd"]],
        direction = params[["direction"]],
        x_dist = params[["x_dist"]],
        x_prevalence = params[["x_prevalence"]],
        x_size = params[["x_size"]],
        x_disp_param = params[["x_disp_param"]],
        rand_intercept_sd = params[["rand_intercept_sd"]],
        rand_slope_sd = params[["rand_slope_sd"]],
        rand_eff_corr = params[["rand_eff_corr"]],
        gamma0 = params[["gamma0"]],
        gamma1 = params[["gamma1"]],
        gamma2 = params[["gamma2"]],
        gamma_sd = params[["gamma_sd"]],
        time_scale = params[["time_scale"]]
    )
    return(data)
}
