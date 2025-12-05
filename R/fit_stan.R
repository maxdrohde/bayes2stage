#' Fit a Bayesian two-stage model using Stan
#'
#' Fits a mixed effects model with imputation using Stan via the instantiate package.
#'
#' @param data A data frame containing the outcome and covariates
#' @param main_model_covariates Character vector of covariate names for the main model
#' @param imputation_model_covariates Character vector of covariate names for the imputation model
#' @param nchains Number of MCMC chains (default: 4)
#' @param iter_warmup Number of warmup iterations per chain (default: 1000)
#' @param iter_sampling Number of sampling iterations per chain (default: 1000)
#' @param adapt_delta Target acceptance rate for HMC (default: 0.8)
#' @param seed Random seed for reproducibility (default: 777L)
#' @param parallel_chains Number of chains to run in parallel (default: 1L)
#' @return A CmdStanMCMC fit object
#' @export
fit_stan_model <- function(data,
                           main_model_covariates,
                           imputation_model_covariates,
                           nchains = 4,
                           iter_warmup = 1000,
                           iter_sampling = 1000,
                           adapt_delta = 0.8,
                           seed = 777L,
                           parallel_chains = 1L) {

  data_list <- format_data_mcmc(data,
                                main_vars = main_model_covariates,
                                imputation_vars = imputation_model_covariates)

  mod <- instantiate::stan_package_model(
    name = "mixed_effects_imputation",
    package = "bayes2stage"
  )

  fit <- mod$sample(
    data = data_list,
    seed = seed,
    chains = nchains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta
  )

  return(fit)
}
