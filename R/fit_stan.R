#' Fit a Bayesian two-stage model using Stan
#'
#' Fits a mixed effects model with imputation using Stan via the instantiate package.
#'
#' @param data A data frame containing the outcome and covariates
#' @param main_model_covariates Character vector of covariate names for the main model
#' @param imputation_model_covariates Character vector of covariate names for the imputation model
#' @param imputation_distribution Distribution for the imputation model:
#'   "normal" for continuous x, "bernoulli" for binary x, "beta_binomial"
#'   for bounded count data, or "negative_binomial" for unbounded count data
#'   (default: "normal")
#' @param n_chains Number of MCMC chains (default: 4)
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
                           imputation_distribution = c("normal",
                                                       "bernoulli",
                                                       "beta_binomial",
                                                       "negative_binomial"),
                           n_chains = 4,
                           iter_warmup = 1000,
                           iter_sampling = 1000,
                           adapt_delta = 0.8,
                           seed = 777L,
                           parallel_chains = 1L) {

  imputation_distribution <- match.arg(imputation_distribution)

  data_list <- format_data_mcmc(data,
                                main_model_covariates = main_model_covariates,
                                imputation_model_covariates = imputation_model_covariates,
                                imputation_distribution = imputation_distribution)

  model_name <- paste0("mixed_effects_imputation_", imputation_distribution)

  mod <- instantiate::stan_package_model(
    name = model_name,
    package = "bayes2stage"
  )

  fit <- mod$sample(
    data = data_list,
    seed = seed,
    chains = n_chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta
  )

  return(fit)
}
