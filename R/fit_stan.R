#' Run Stan inference with error handling
#'
#' @param method_fn Function to call for inference
#' @param method_name Name of the method for error messages
#' @param start_msg Message to display when starting
#' @param success_msg Message to display on success
#' @param warning_msg Optional warning message to display after success
#' @return The fit object from the inference method
#' @noRd
run_stan_inference <- function(method_fn,
                               method_name,
                               start_msg,
                               success_msg,
                               warning_msg = NULL) {
    cli::cli_alert_info(start_msg)
    fit <- tryCatch(
        method_fn(),
        error = function(e) {
            cli::cli_abort(
                c("{method_name} failed:",
                  "x" = e$message),
                call = NULL
            )
        }
    )
    cli::cli_alert_success(success_msg)
    if (!is.null(warning_msg)) {
        cli::cli_alert_warning(warning_msg)
    }
    return(fit)
}

#' Execute Pathfinder algorithm
#'
#' Shared helper for running Pathfinder, used both as standalone method
#' and for MCMC initialization.
#'
#' @param mod CmdStan model object
#' @param data_list Data list for Stan
#' @param seed Random seed
#' @param num_paths Number of Pathfinder paths
#' @param draws Number of draws (NULL for initialization use)
#' @return CmdStanPathfinder object
#' @noRd
execute_pathfinder <- function(mod, data_list, seed, num_paths, draws = NULL) {
    args <- list(
        data = data_list,
        seed = seed,
        num_paths = num_paths
    )
    if (!is.null(draws)) {
        args$draws <- draws
    }
    result <- do.call(mod$pathfinder, args)
    return(result)
}

#' Run Pathfinder variational inference
#' @noRd
run_pathfinder <- function(mod, data_list, seed, num_paths, draws) {
    run_stan_inference(
        method_fn = function() {
            execute_pathfinder(mod, data_list, seed, num_paths, draws)
        },
        method_name = "Pathfinder",
        start_msg = "Running Pathfinder variational inference...",
        success_msg = "Pathfinder complete."
    )
}

#' Run Laplace approximation
#' @noRd
run_laplace <- function(mod, data_list, seed, draws) {
    run_stan_inference(
        method_fn = function() {
            mod$laplace(
                data = data_list,
                seed = seed,
                draws = draws
            )
        },
        method_name = "Laplace approximation",
        start_msg = "Running Laplace approximation...",
        success_msg = "Laplace approximation complete."
    )
}

#' Run optimization (MAP estimate)
#' @noRd
run_optimize <- function(mod, data_list, seed, algorithm, iter) {
    run_stan_inference(
        method_fn = function() {
            mod$optimize(
                data = data_list,
                seed = seed,
                algorithm = algorithm,
                iter = iter
            )
        },
        method_name = "Optimization",
        start_msg = "Running optimization (MAP estimate)...",
        success_msg = "Optimization complete.",
        warning_msg = "Note: Optimize returns point estimates only, not posterior samples."
    )
}

#' Run MCMC sampling
#' @noRd
run_mcmc <- function(mod, data_list, seed, n_chains, parallel_chains,
                     iter_warmup, iter_sampling, adapt_delta,
                     use_pathfinder_init, pathfinder_num_paths) {
    init <- NULL
    if (use_pathfinder_init) {
        init <- tryCatch(
            {
                cli::cli_alert_info("Running Pathfinder for initialization...")
                pathfinder_fit <- execute_pathfinder(
                    mod, data_list, seed, pathfinder_num_paths
                )
                cli::cli_alert_success("Pathfinder complete. Initializing MCMC chains.")
                pathfinder_fit
            },
            error = function(e) {
                cli::cli_warn(c(
                    "Pathfinder initialization failed.",
                    "x" = e$message,
                    "i" = "Using default initialization instead."
                ))
                return(NULL)
            }
        )
    }

    run_stan_inference(
        method_fn = function() {
            mod$sample(
                data = data_list,
                seed = seed,
                init = init,
                chains = n_chains,
                parallel_chains = parallel_chains,
                iter_warmup = iter_warmup,
                iter_sampling = iter_sampling,
                adapt_delta = adapt_delta
            )
        },
        method_name = "MCMC sampling",
        start_msg = "Running MCMC sampling...",
        success_msg = "MCMC sampling complete."
    )
}

#' Fit a Bayesian two-stage model using Stan
#'
#' Fits a mixed effects model with imputation using Stan via the instantiate package.
#'
#' @param data A data frame containing the outcome and covariates
#' @param main_model_formula One-sided formula or string for covariates in the
#'   main model (e.g., `~ age + splines::ns(bmi, 3)`). Intercept is
#'   automatically removed.
#' @param imputation_model_formula One-sided formula or string for covariates
#'   in the imputation model (e.g., `~ age + factor(site)`). Intercept is
#'   automatically removed.
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
#'   - "marginalized_mixture": Like "marginalized" but models x with a mixture of
#'     normals for multimodal distributions. Only available for normal imputation.
#'     Requires specifying `mixture_components`.
#' @param inference_method Inference method to use:
#'   - "mcmc" (default): Full MCMC sampling via NUTS. Most accurate but slowest.
#'   - "pathfinder": Variational inference via Pathfinder algorithm. Fast approximate
#'     inference with posterior draws. Good for initial exploration.
#'   - "laplace": Laplace approximation. Finds posterior mode and approximates with
#'     Gaussian. Very fast with good accuracy for well-behaved posteriors.
#'   - "optimize": Maximum a posteriori (MAP) point estimate only. Fastest but no
#'     uncertainty quantification. Returns CmdStanMLE object (no posterior draws).
#' @param pathfinder_draws Number of approximate posterior draws from Pathfinder
#'   (default: 1000L)
#' @param pathfinder_num_paths Number of Pathfinder paths to run (default: 4L)
#' @param laplace_draws Number of approximate posterior draws from Laplace
#'   (default: 1000L)
#' @param optimize_algorithm Optimization algorithm for MAP estimation:
#'   "lbfgs" (default), "bfgs", or "newton"
#' @param optimize_iter Maximum iterations for optimization (default: 2000L)
#' @param mixture_components Number of mixture components for the imputation model
#'   when using `parameterization = "marginalized_mixture"` (default: 3L)
#' @param use_pathfinder_init Logical; if TRUE, use Pathfinder variational inference
#'   to initialize MCMC chains. Only applies when `inference_method = "mcmc"`.
#'   This can dramatically improve sampling efficiency for complex models with
#'   many latent parameters (default: FALSE)
#' @param n_chains Number of MCMC chains (default: 4)
#' @param iter_warmup Number of warmup iterations per chain (default: 1000)
#' @param iter_sampling Number of sampling iterations per chain (default: 1000)
#' @param adapt_delta Target acceptance rate for HMC (default: 0.8)
#' @param seed Random seed for reproducibility (default: 777L)
#' @param parallel_chains Number of chains to run in parallel (default: 1L)
#' @return A CmdStan fit object. Type depends on `inference_method`:
#'   CmdStanMCMC (mcmc), CmdStanPathfinder (pathfinder), CmdStanLaplace (laplace),
#'   or CmdStanMLE (optimize). All support `$summary()`, but only mcmc/pathfinder/laplace
#'   support `$draws()` for posterior samples.
#' @export
fit_stan_model <- function(data,
                           main_model_formula,
                           imputation_model_formula,
                           imputation_distribution = c("normal",
                                                       "bernoulli",
                                                       "beta_binomial",
                                                       "negative_binomial"),
                           parameterization = c("noncentered",
                                                "centered",
                                                "marginalized",
                                                "marginalized_mixture"),
                           inference_method = c("mcmc",
                                                "pathfinder",
                                                "laplace",
                                                "optimize"),
                           pathfinder_draws = 1000L,
                           pathfinder_num_paths = 4L,
                           laplace_draws = 1000L,
                           optimize_algorithm = c("lbfgs", "bfgs", "newton"),
                           optimize_iter = 2000L,
                           mixture_components = 3L,
                           use_pathfinder_init = FALSE,
                           n_chains = 4L,
                           iter_warmup = 1000L,
                           iter_sampling = 1000L,
                           adapt_delta = 0.8,
                           seed = 777L,
                           parallel_chains = 1L) {

    imputation_distribution <- match.arg(imputation_distribution)
    parameterization <- match.arg(parameterization)
    inference_method <- match.arg(inference_method)
    optimize_algorithm <- match.arg(optimize_algorithm)

    data_list <- format_data_mcmc(data,
                                  main_model_formula = main_model_formula,
                                  imputation_model_formula = imputation_model_formula,
                                  imputation_distribution = imputation_distribution)

    # Validate marginalized_mixture requirements
    if (parameterization == "marginalized_mixture") {
        if (imputation_distribution != "normal") {
            cli::cli_abort(c(
                "`parameterization = \"marginalized_mixture\"` only supports normal imputation.",
                "x" = "Current imputation_distribution: {imputation_distribution}",
                "i" = "Use `imputation_distribution = \"normal\"` with marginalized_mixture."
            ))
        }
        if (mixture_components < 2L) {
            cli::cli_abort("`mixture_components` must be at least 2.")
        }
        data_list$K <- as.integer(mixture_components)
    }

    model_name <- glue::glue("mixed_effects_imputation_{imputation_distribution}")
    if (parameterization == "centered") {
        model_name <- glue::glue("{model_name}_centered")
    } else if (parameterization == "marginalized") {
        model_name <- glue::glue("{model_name}_marginalized")
    } else if (parameterization == "marginalized_mixture") {
        model_name <- glue::glue("{model_name}_marginalized_mixture")
    }

    mod <- instantiate::stan_package_model(
        name = model_name,
        package = "bayes2stage"
    )

    if (inference_method != "mcmc" && use_pathfinder_init) {
        cli::cli_warn("`use_pathfinder_init` only applies to `inference_method = \"mcmc\"`")
    }

    fit <- switch(
        inference_method,
        pathfinder = run_pathfinder(mod, data_list, seed,
                                    pathfinder_num_paths, pathfinder_draws),
        laplace = run_laplace(mod, data_list, seed, laplace_draws),
        optimize = run_optimize(mod, data_list, seed, optimize_algorithm, optimize_iter),
        mcmc = run_mcmc(mod, data_list, seed, n_chains, parallel_chains,
                        iter_warmup, iter_sampling, adapt_delta,
                        use_pathfinder_init, pathfinder_num_paths)
    )

    return(fit)
}
