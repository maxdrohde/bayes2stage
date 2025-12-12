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
  stopifnot(i >= 1L, j >= 1L)
  k <- i + j
  seed <- as.integer((k * (k + 1L)) %/% 2L + j)
  if (seed > .Machine$integer.max) stop("seed overflow")
  seed
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

#' Extract MCMC summary statistics
#'
#' A wrapper around MCMCvis::MCMCsummary to extract summary statistics from
#' MCMC samples. For CmdStanR fits, also extracts HMC diagnostics (divergent
#' transitions, max treedepth exceeded, E-BFMI).
#'
#' @param mcmc_output MCMC output object
#' @param dataset_id A character string identifying the dataset (added as a column)
#' @return A data frame with summary statistics for each parameter
#' @export
mcmc_summary <- function(mcmc_output,
                         dataset_id){

  out <-
  mcmc_output |>
    MCMCvis::MCMCsummary() |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "parameter")

  out$dataset_id <- dataset_id


  # Extract HMC diagnostics if this is a CmdStanMCMC object
  if (inherits(mcmc_output, "CmdStanMCMC")) {
    diag <- mcmc_output$diagnostic_summary()
    out$divergent_transitions <- sum(diag$num_divergent)
    out$max_treedepth_exceeded <- sum(diag$num_max_treedepth)
    out$ebfmi_min <- min(diag$ebfmi)
  } else {
    out$divergent_transitions <- NA_integer_
    out$max_treedepth_exceeded <- NA_integer_
    out$ebfmi_min <- NA_real_
  }

  return(out)
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
  if (length(missing)) {
    stop("`data` is missing required column",
         if (length(missing) > 1) "s" else "", ": ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}
