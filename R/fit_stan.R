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
#' @param parameterization Parameterization for random effects:
#'   - "noncentered" (default): works better with weakly informative data (small N)
#'   - "centered": works better with highly informative data (large N)
#'   - "marginal": integrates out random effects analytically - RECOMMENDED for
#'     large N (>500) as it dramatically reduces the number of parameters and
#'     improves mixing. Only available for normal imputation distribution.
#'   - "improved": centered parameterization with better-scaled priors
#' @param n_chains Number of MCMC chains (default: 4)
#' @param iter_warmup Number of warmup iterations per chain (default: 1000)
#' @param iter_sampling Number of sampling iterations per chain (default: 1000)
#' @param adapt_delta Target acceptance rate for HMC (default: 0.8). For large N
#'   or complex posteriors, consider increasing to 0.95 or 0.99.
#' @param max_treedepth Maximum tree depth for NUTS (default: 10). Increase to
#'   12-15 if you see max treedepth warnings.
#' @param seed Random seed for reproducibility (default: 777L)
#' @param parallel_chains Number of chains to run in parallel (default: 1L)
#' @param prior_beta_sd SD for normal priors on regression coefficients (default: 5).
#'   Only used for "marginal" and "improved" parameterizations.
#' @param prior_sigma_rate Rate for exponential prior on residual SD (default: 1).
#'   Only used for "marginal" and "improved" parameterizations.
#' @param prior_sigma_re_rate Rate for exponential prior on random effect SDs
#'   (default: 0.5). Only used for "marginal" and "improved" parameterizations.
#' @return A CmdStanMCMC fit object
#' @export
fit_stan_model <- function(data,
                           main_model_covariates,
                           imputation_model_covariates,
                           imputation_distribution = c("normal",
                                                       "bernoulli",
                                                       "beta_binomial",
                                                       "negative_binomial"),
                           parameterization = c("noncentered", "centered",
                                               "marginal", "improved"),
                           n_chains = 4,
                           iter_warmup = 1000,
                           iter_sampling = 1000,
                           adapt_delta = 0.8,
                           max_treedepth = 10,
                           seed = 777L,
                           parallel_chains = 1L,
                           prior_beta_sd = 5.0,
                           prior_sigma_rate = 1.0,
                           prior_sigma_re_rate = 0.5) {

  imputation_distribution <- match.arg(imputation_distribution)
  parameterization <- match.arg(parameterization)

  # Validate parameterization choices

  if (parameterization %in% c("marginal", "improved") &&
      imputation_distribution != "normal") {
    stop("The '", parameterization, "' parameterization is only available for ",
         "normal imputation distribution. Use 'centered' or 'noncentered' for ",
         "other distributions.")
  }

  # Format data based on parameterization
  data_list <- format_data_mcmc(data,
                                main_model_covariates = main_model_covariates,
                                imputation_model_covariates = imputation_model_covariates,
                                imputation_distribution = imputation_distribution)

  # Add prior hyperparameters for marginal and improved models
  if (parameterization %in% c("marginal", "improved")) {
    data_list$prior_beta_sd <- prior_beta_sd
    data_list$prior_sigma_rate <- prior_sigma_rate
    data_list$prior_sigma_re_rate <- prior_sigma_re_rate
  }

  # Determine model name
  model_name <- paste0("mixed_effects_imputation_", imputation_distribution)
  if (parameterization == "centered") {
    model_name <- paste0(model_name, "_centered")
  } else if (parameterization == "marginal") {
    model_name <- paste0(model_name, "_marginal")
  } else if (parameterization == "improved") {
    model_name <- paste0(model_name, "_improved")
  }

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
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
  )

  return(fit)
}
