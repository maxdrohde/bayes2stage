is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

is_positive_integer <- function(x) {
  is.numeric(x) & is.finite(x) & x > 0 & is.wholenumber(x)
}

#' @export
mcmc_forest <- function(mcmc_output){
  MCMCvis::MCMCplot(mcmc_output)

  if (print_to_pdf) {
    MCMCvis::MCMCtrace(mcmc_output)
  }
}

#' Plot / Summarize MCMC Output
#' @export
mcmc_trace <- function(mcmc_output,
                       print_to_pdf = FALSE){

    MCMCvis::MCMCtrace(mcmc_output,
                       pdf = print_to_pdf)

}

#' Extract MCMC summary
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

#' @export
check_cols <- function(df, required_cols) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing)) {
    stop("`df` is missing required column",
         if (length(missing) > 1) "s" else "", ": ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}
