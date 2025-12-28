is_whole_number <- function(x, tol = .Machine$double.eps^0.5) {
    result <- abs(x - round(x)) < tol
    return(result)
}

is_positive_integer <- function(x) {
    result <- is.numeric(x) & is.finite(x) & x > 0 & is_whole_number(x)
    return(result)
}

#' Generate a unique seed from two indices using Cantor pairing
#'
#' Uses the Cantor pairing function to generate a unique integer seed from
#' two positive integer indices. Useful for reproducible simulations where
#' each (i, j) combination needs a unique seed.
#'
#' @param i First positive integer index
#' @param j Second positive integer index
#' @return A unique integer seed
#' @examples
#' cantor_seed(1, 1)
#' cantor_seed(5, 10)
#' @export
cantor_seed <- function(i, j) {
  if (i < 1L || j < 1L) {
    cli::cli_abort("{.arg i} and {.arg j} must be >= 1.")
  }
  k <- i + j
  seed <- as.integer((k * (k + 1L)) %/% 2L + j)
  if (seed > .Machine$integer.max) {
    cli::cli_abort("Seed overflow: generated seed exceeds integer max.")
  }
  return(seed)
}

#' Create a forest plot of MCMC output
#'
#' A wrapper around MCMCvis::MCMCplot to create forest plots of MCMC samples.
#'
#' @param mcmc_output MCMC output object
#' @return A forest plot of parameter estimates
#' @export
mcmc_forest <- function(mcmc_output) {
    p <- MCMCvis::MCMCplot(mcmc_output)
    return(p)
}

#' Plot MCMC trace plots
#'
#' A wrapper around MCMCvis::MCMCtrace to create trace plots of MCMC samples.
#'
#' @param mcmc_output MCMC output object
#' @param print_to_pdf Logical; if TRUE, saves trace plots to a PDF file
#' @return Trace plots of MCMC samples
#' @export
mcmc_trace <- function(mcmc_output, print_to_pdf = FALSE) {
    p <- MCMCvis::MCMCtrace(mcmc_output, pdf = print_to_pdf)
    return(p)
}

#' Extract summary statistics from model fit objects
#'
#' S3 generic for extracting summary statistics from various fit object types.
#' Supports CmdStanR fit types and ACML fits.
#'
#' @param object A model fit object
#' @param ... Additional arguments passed to methods
#' @return A data frame with summary statistics for each parameter
#' @export
model_summary <- function(object, ...) {
    UseMethod("model_summary")
}

#' @rdname model_summary
#' @export
model_summary.CmdStanMCMC <- function(object, ...) {
    out <- object$summary(
        variables = get_key_parameters(),
        mean, median, sd,
        ~quantile(.x, probs = c(0.025, 0.05, 0.5, 0.95, 0.975)),
        posterior::rhat, posterior::ess_bulk, posterior::ess_tail
    ) |>
        as.data.frame() |>
        dplyr::rename(
            parameter = variable,
            q2_5 = `2.5%`,
            q5 = `5%`,
            q50 = `50%`,
            q95 = `95%`,
            q97_5 = `97.5%`,
            rhat = `posterior::rhat`,
            ess_bulk = `posterior::ess_bulk`,
            ess_tail = `posterior::ess_tail`
        )

    diag <- object$diagnostic_summary()
    out$divergent_transitions <- sum(diag$num_divergent, na.rm = TRUE)
    out$max_treedepth_exceeded <- sum(diag$num_max_treedepth, na.rm = TRUE)
    ebfmi_vals <- diag$ebfmi[is.finite(diag$ebfmi)]
    out$ebfmi_min <- if (length(ebfmi_vals) > 0L) min(ebfmi_vals) else NA_real_

    return(out)
}

#' @rdname model_summary
#' @export
model_summary.CmdStanPathfinder <- function(object, ...) {
    out <- object$summary(
        variables = get_key_parameters(),
        mean, median, sd,
        ~quantile(.x, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))
    ) |>
        as.data.frame() |>
        dplyr::rename(
            parameter = variable,
            q2_5 = `2.5%`,
            q5 = `5%`,
            q50 = `50%`,
            q95 = `95%`,
            q97_5 = `97.5%`
        )

    # Pathfinder doesn't have MCMC diagnostics
    out$rhat <- NA_real_
    out$ess_bulk <- NA_real_
    out$ess_tail <- NA_real_
    out$divergent_transitions <- NA_integer_
    out$max_treedepth_exceeded <- NA_integer_
    out$ebfmi_min <- NA_real_

    return(out)
}

#' @rdname model_summary
#' @export
model_summary.CmdStanLaplace <- function(object, ...) {
    out <- object$summary(
        variables = get_key_parameters(),
        mean, median, sd,
        ~quantile(.x, probs = c(0.025, 0.05, 0.5, 0.95, 0.975))
    ) |>
        as.data.frame() |>
        dplyr::rename(
            parameter = variable,
            q2_5 = `2.5%`,
            q5 = `5%`,
            q50 = `50%`,
            q95 = `95%`,
            q97_5 = `97.5%`
        )

    # Laplace doesn't have MCMC diagnostics
    out$rhat <- NA_real_
    out$ess_bulk <- NA_real_
    out$ess_tail <- NA_real_
    out$divergent_transitions <- NA_integer_
    out$max_treedepth_exceeded <- NA_integer_
    out$ebfmi_min <- NA_real_

    return(out)
}

#' @rdname model_summary
#' @export
model_summary.CmdStanMLE <- function(object, ...) {
    summ <- object$summary() |>
        as.data.frame()

    key_pars <- get_key_parameters()
    summ <- summ[summ$variable %in% key_pars, ]

    out <- data.frame(
        parameter = summ$variable,
        mean = summ$estimate,
        median = summ$estimate,
        sd = NA_real_,
        q2_5 = NA_real_,
        q5 = NA_real_,
        q50 = summ$estimate,
        q95 = NA_real_,
        q97_5 = NA_real_,
        rhat = NA_real_,
        ess_bulk = NA_real_,
        ess_tail = NA_real_,
        check.names = FALSE
    )
    out$divergent_transitions <- NA_integer_
    out$max_treedepth_exceeded <- NA_integer_
    out$ebfmi_min <- NA_real_

    return(out)
}

#' @rdname model_summary
#' @export
model_summary.acml_fit <- function(object, ...) {
    out <- object
    if ("variable" %in% names(out) && !"parameter" %in% names(out)) {
        out <- dplyr::rename(out, parameter = variable)
    }
    out$parameter <- translate_parameter_names(out$parameter, from = "acml", to = "stan")

    # Compute frequentist 95% CI from normal approximation
    z_975 <- qnorm(0.975)
    out$q2_5 <- out$mean - z_975 * out$sd
    out$q97_5 <- out$mean + z_975 * out$sd

    # Other quantiles from normal approximation
    out$q5 <- out$mean - qnorm(0.95) * out$sd
    out$q50 <- out$mean
    out$q95 <- out$mean + qnorm(0.95) * out$sd

    # MCMC diagnostics not applicable
    out$rhat <- NA_real_
    out$ess_bulk <- NA_real_
    out$ess_tail <- NA_real_
    out$divergent_transitions <- NA_integer_
    out$max_treedepth_exceeded <- NA_integer_
    out$ebfmi_min <- NA_real_

    return(out)
}

#' @rdname model_summary
#' @export
model_summary.default <- function(object, ...) {
    cli::cli_abort("No {.fn model_summary} method for class {.cls {class(object)}}")
}

#' Check for required columns in a data frame
#'
#' Validates that a data frame contains all required columns.
#'
#' @param data A data frame to check
#' @param required_cols A character vector of required column names
#' @return Invisibly returns TRUE if all columns are present; otherwise throws an error
#' @export
check_cols <- function(data, required_cols) {
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0L) {
    cli::cli_abort("{.arg data} is missing required column{?s}: {.val {missing}}.")
  }
  invisible(TRUE)
}

validate_sampling_type <- function(sampling_type) {
  if (!(sampling_type %in% c("intercept", "slope"))) {
    cli::cli_abort("{.arg sampling_type} must be {.val intercept} or {.val slope}.")
  }
  invisible(TRUE)
}
