#' Prepare a formula for model matrix creation
#'
#' Converts string to formula if needed and validates input. The intercept is
#' kept in the formula so that `model.matrix()` uses K-1 (treatment) coding for
#' factor variables. The intercept column is dropped after `model.matrix()` in
#' [format_data_mcmc()] because the Stan model has its own intercept parameter.
#'
#' @param x A formula or character string
#' @param name Name of the argument (for error messages)
#' @return A validated formula (intercept retained for proper factor encoding)
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
    x
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

  # Validate required columns

  required_cols <- c("id", "t", "y", "x")
  if (imputation_distribution == "beta_binomial") {
      required_cols <- c(required_cols, "n_trials")
  }
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0L) {
      cli::cli_abort("Missing required columns in {.arg data}: {.field {missing_cols}}")
  }

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

  stopifnot(
      "`main_model_formula` must not reference `x` or `t` (modeled as separate Stan parameters)." =
          !any(c("x", "t") %in% all.vars(main_model_formula))
  )

  X <- stats::model.matrix(main_model_formula, data = data, na.action = stats::na.pass)
  Z <- stats::model.matrix(imputation_model_formula, data = id_df, na.action = stats::na.pass)

  # Drop the intercept column — the Stan model has its own alpha_main / alpha_imputation.
  # We build model.matrix WITH the intercept so factor variables get proper K-1
  # (treatment) encoding rather than K indicator columns.
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  Z <- Z[, colnames(Z) != "(Intercept)", drop = FALSE]

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
