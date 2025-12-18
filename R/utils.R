is_whole_number <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

is_positive_integer <- function(x) {
  is.numeric(x) & is.finite(x) & x > 0 & is_whole_number(x)
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
mcmc_forest <- function(mcmc_output){
  MCMCvis::MCMCplot(mcmc_output)
}

#' Plot MCMC trace plots
#'
#' A wrapper around MCMCvis::MCMCtrace to create trace plots of MCMC samples.
#'
#' @param mcmc_output MCMC output object
#' @param print_to_pdf Logical; if TRUE, saves trace plots to a PDF file
#' @return Trace plots of MCMC samples
#' @export
mcmc_trace <- function(mcmc_output,
                       print_to_pdf = FALSE){

    MCMCvis::MCMCtrace(mcmc_output,
                       pdf = print_to_pdf)

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
    out <- object$summary() |>
        as.data.frame() |>
        dplyr::rename(parameter = "variable")

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
    out <- object$summary() |>
        as.data.frame() |>
        dplyr::rename(parameter = "variable")

    out$divergent_transitions <- NA_integer_
    out$max_treedepth_exceeded <- NA_integer_
    out$ebfmi_min <- NA_real_

    return(out)
}

#' @rdname model_summary
#' @export
model_summary.CmdStanLaplace <- function(object, ...) {
    out <- object$summary() |>
        as.data.frame() |>
        dplyr::rename(parameter = "variable")

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

    out <- data.frame(
        parameter = summ$variable,
        mean = summ$estimate,
        sd = NA_real_,
        median = summ$estimate,
        q5 = NA_real_,
        q95 = NA_real_,
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
    # Rename variable to parameter if needed
    if ("variable" %in% names(out) && !"parameter" %in% names(out)) {
        out <- dplyr::rename(out, parameter = "variable")
    }
    out$parameter <- translate_parameter_names(out$parameter, from = "acml", to = "stan")
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
