is_whole_number <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

is_positive_integer <- function(x) {
  is.numeric(x) & is.finite(x) & x > 0 & is_whole_number(x)
}

#' Create a forest plot of MCMC output
#'
#' A wrapper around MCMCvis::MCMCplot to create forest plots of MCMC samples.
#'
#' @param mcmc_output MCMC output object (e.g., from fit_model)
#' @return A forest plot of parameter estimates
#' @export
mcmc_forest <- function(mcmc_output){
  MCMCvis::MCMCplot(mcmc_output)
}

#' Plot MCMC trace plots
#'
#' A wrapper around MCMCvis::MCMCtrace to create trace plots of MCMC samples.
#'
#' @param mcmc_output MCMC output object (e.g., from fit_model)
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
#' MCMC samples.
#'
#' @param mcmc_output MCMC output object (e.g., from fit_model)
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
