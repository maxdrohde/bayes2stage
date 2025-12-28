#' Prepare a formula for model matrix creation
#'
#' Converts string to formula if needed, validates input, and removes intercept.
#'
#' @param x A formula or character string
#' @param name Name of the argument (for error messages)
#' @return A formula with intercept suppressed
#' @noRd
prepare_formula <- function(x, name) {
    if (is.null(x)) {
        cli::cli_abort("{.arg {name}} must be provided")
    }
    if (is.character(x)) {
        x <- stats::as.formula(x)
    }
    if (!inherits(x, "formula")) {
        cli::cli_abort("{.arg {name}} must be a formula or string, not {.cls {class(x)}}")
    }
    stats::update(x, ~ . - 1)
}

#' Format the simulated data for Stan
#'
#' @param data Dataset to use
#' @param main_model_formula One-sided formula or string for covariates in the
#'   main model (e.g., `~ age + splines::ns(bmi, 3)`). Intercept is
#'   automatically removed.
#' @param imputation_model_formula One-sided formula or string for covariates
#'   in the imputation model (e.g., `~ age + factor(site)`). Intercept is
#'   automatically removed.
#' @param imputation_distribution Distribution for the imputation model:
#'   "normal" for continuous x, "bernoulli" for binary x, "beta_binomial"
#'   for bounded count data, or "negative_binomial" for unbounded count data
#' @return A list suitable for input to MCMC software
#' @export
format_data_mcmc <- function(data,
                             main_model_formula = NULL,
                             imputation_model_formula = NULL,
                             imputation_distribution = c("normal",
                                                         "bernoulli",
                                                         "beta_binomial",
                                                         "negative_binomial"))
  {

  imputation_distribution <- match.arg(imputation_distribution)

  data <-
    data |>
    dplyr::arrange(id, t)

  # Re-index id to 1..G
  data <-
    data |>
    dplyr::mutate(id_idx = as.integer(factor(id)))

  id_df <-
    data |>
    dplyr::distinct(id_idx,
                    .keep_all = TRUE) |>
    dplyr::arrange(id_idx)

  main_model_formula <- prepare_formula(main_model_formula, "main_model_formula")
  imputation_model_formula <- prepare_formula(imputation_model_formula, "imputation_model_formula")

  X <- stats::model.matrix(main_model_formula, data = data, na.action = stats::na.pass)
  Z <- stats::model.matrix(imputation_model_formula, data = id_df, na.action = stats::na.pass)

  P <- ncol(X)
  S <- ncol(Z)

  # Format x_obs based on imputation_distribution
  x_obs_raw <- id_df$x[!is.na(id_df$x)]
  if (imputation_distribution %in% c("bernoulli",
                                     "beta_binomial",
                                     "negative_binomial")) {
    x_obs <- as.integer(x_obs_raw)
  } else {
    x_obs <- x_obs_raw
  }

  # Compute pos and len arrays for optimized Stan models
  # These map from subject IDs (G) to observation indices (N)
  G <- nrow(id_df)
  pos <- integer(G)
  len <- integer(G)
  for (g in seq_len(G)) {
    idx <- which(data$id_idx == g)
    pos[g] <- min(idx)
    len[g] <- length(idx)
  }

  data_list <-
    list(
      N       = nrow(data),
      G       = G,
      G_obs   = sum(!is.na(id_df$x)),
      G_mis   = sum(is.na(id_df$x)),
      P       = P,
      S       = S,
      t       = data$t,
      X       = X,
      Z       = Z,
      y       = data$y,
      index_obs = which(!is.na(id_df$x)),
      index_mis = which(is.na(id_df$x)),
      x_obs     = x_obs,
      id        = data$id_idx,
      pos       = pos,
      len       = len
    )

  # Add n_trials for beta_binomial model (array, same interface for standard and marginalized)
  if (imputation_distribution == "beta_binomial") {
    data_list$n_trials <- as.integer(id_df$n_trials)
  }

  # Add x_max for negative_binomial model (same name for standard and marginalized)
  if (imputation_distribution == "negative_binomial") {
    # Use 3x the max observed value or 100, whichever is larger
    max_observed <- if (length(x_obs) > 0) max(x_obs) else 0
    data_list$x_max <- as.integer(max(max_observed * 3, 100))
  }

  return(data_list)
}
