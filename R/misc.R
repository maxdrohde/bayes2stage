is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

#' Plot / Summarize MCMC Output
#' @export
check_mcmc <- function(mcmc_output){
  print(MCMCvis::MCMCsummary(mcmc_output))
  MCMCvis::MCMCtrace(mcmc_output)
  MCMCvis::MCMCplot(mcmc_output)
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
