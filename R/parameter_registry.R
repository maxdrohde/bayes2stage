#' @title Parameter Registry
#' @name parameter_registry
#' @description Central registry for parameter name mappings across contexts
NULL

#' Get parameter registry with all name mappings
#'
#' Central registry mapping parameter names across different contexts:
#' grid names (simulation config), Stan output names, display names, and ACML names.
#'
#' @return A tibble with columns: grid_name, stan_name, display_name, acml_name
#' @export
get_parameter_registry <- function() {
    # Note: ACML reports SDs on log scale and correlation on Fisher z scale,
    # so variance components are not directly comparable to Stan output
    result <- tibble::tribble(
        ~grid_name,             ~stan_name,             ~display_name,            ~acml_name,
        "beta_x",               "beta_x",               "Beta x",                 "beta_x",
        "beta_z",               "beta[1]",              "Beta z",                 "beta_z",
        "beta_t",               "beta_t",               "Beta t",                 "beta_t",
        "beta_t_x_interaction", "beta_t_x_interaction", "Beta x:t Interaction",   "beta_t_x_interaction",
        "alpha_main",           "alpha_main",           "Intercept (Main)",       "alpha_main",
        "error_sd",             "sigma_main",           "Error SD (Main)",        NA_character_,
        "rand_intercept_sd",    "sigma_re[1]",          "Random Intercept SD",    NA_character_,
        "rand_slope_sd",        "sigma_re[2]",          "Random Slope SD",        NA_character_,
        "gamma0",               "alpha_imputation",     "Intercept (Imputation)", NA_character_,
        "gamma1",               "gamma[1]",             "Gamma_1",                NA_character_,
        "gamma2",               "gamma[2]",             "Gamma_2",                NA_character_,
        "gamma_sd",             "sigma_imputation",     "Error SD (Imputation)",  NA_character_
    )
    return(result)
}

#' Translate parameter names between contexts
#'
#' @param names Character vector of parameter names to translate
#' @param from Source context: "grid", "stan", "display", or "acml"
#' @param to Target context: "grid", "stan", "display", or "acml"
#' @return Character vector of translated names (preserves original if no mapping found)
#' @export
translate_parameter_names <- function(names, from = "stan", to = "display") {
    registry <- get_parameter_registry()
    from_col <- glue::glue("{from}_name")
    to_col <- glue::glue("{to}_name")

    valid_cols <- c("grid_name", "stan_name", "display_name", "acml_name")
    if (!from_col %in% valid_cols) {
        cli::cli_abort("{.arg from} must be one of: grid, stan, display, acml")
    }
    if (!to_col %in% valid_cols) {
        cli::cli_abort("{.arg to} must be one of: grid, stan, display, acml")
    }

    idx <- match(names, registry[[from_col]])
    result <- registry[[to_col]][idx]
    # Keep original name if no translation found
    result <- dplyr::coalesce(result, names)
    return(result)
}

#' Get default key parameters for diagnostics
#'
#' Returns the Stan parameter names that are typically of interest for
#' MCMC diagnostics and model summaries.
#'
#' @return Character vector of Stan parameter names
#' @export
get_key_parameters <- function() {
    registry <- get_parameter_registry()
    result <- registry$stan_name[!is.na(registry$stan_name)]
    return(result)
}
