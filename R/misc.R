is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

#' Plot / Summarize MCMC Output
#' @export
check_mcmc <- function(mcmc_output){
  print(MCMCvis::MCMCsummary(mcmc_output))
  MCMCvis::MCMCtrace(mcmc_output)
  MCMCvis::MCMCplot(mcmc_output)
}
