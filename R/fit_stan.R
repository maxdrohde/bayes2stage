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
#'   - "noncentered" (default): Works better with weakly informative data (small N)
#'   - "centered": Works better with highly informative data (large N)
#'   - "marginalized": Integrates out random effects and missing x analytically.
#'     Eliminates funnel geometry but has O(G * n_g^3) likelihood cost per subject.
#'     Best for datasets with sampling issues. Available for all distributions.
#'     Note: Beta-binomial marginalized requires all subjects have the same n_trials.
#' @param use_pathfinder_init Logical; if TRUE, use Pathfinder variational inference
#'   to initialize MCMC chains. This can dramatically improve sampling efficiency
#'   for complex models with many latent parameters (default: FALSE)
#' @param pathfinder_num_paths Number of Pathfinder paths to run (default: 4, one per chain)
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
                           parameterization = c("noncentered",
                                                "centered",
                                                "marginalized"),
                           use_pathfinder_init = FALSE,
                           pathfinder_num_paths,
                           n_chains,
                           iter_warmup,
                           iter_sampling,
                           adapt_delta = 0.8,
                           seed = 777L,
                           parallel_chains = 1L) {

  imputation_distribution <- match.arg(imputation_distribution)
  parameterization <- match.arg(parameterization)

  data_list <- format_data_mcmc(data,
                                main_model_covariates = main_model_covariates,
                                imputation_model_covariates = imputation_model_covariates,
                                imputation_distribution = imputation_distribution)


  model_name <- glue::glue("mixed_effects_imputation_{imputation_distribution}")
  if (parameterization == "centered") {
    model_name <- glue::glue("{model_name}_centered")
  } else if (parameterization == "marginalized") {
    model_name <- glue::glue("{model_name}_marginalized")
  }

  mod <- instantiate::stan_package_model(
    name = model_name,
    package = "bayes2stage"
  )

  init <- NULL
  if (use_pathfinder_init) {
    cli::cli_alert_info("Running Pathfinder for initialization...")
    pathfinder_fit <- mod$pathfinder(
      data = data_list,
      seed = seed,
      num_paths = pathfinder_num_paths,
      num_threads = parallel_chains
    )
    init <- pathfinder_fit
    cli::cli_alert_success("Pathfinder complete. Initializing MCMC chains.")
  }

  fit <- mod$sample(
    data = data_list,
    seed = seed,
    init = init,
    chains = n_chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    adapt_delta = adapt_delta
  )

  return(fit)
}
