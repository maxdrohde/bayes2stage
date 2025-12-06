#' @import nimble
#' @import nimbleEcology
#' @import nimbleMacros

use_waic <- TRUE

nimbleOptions("enableWAIC" = use_waic)

###############################################################################
# Create the NIMBLE model line for the main model
# Allows the user to input any number of covariates
create_main_model_linpred <- function(main_model_covariates,
                                      main_model_priors){

  expanded_covariates <-
    purrr::map_chr(main_model_covariates, \(x) glue::glue("{x}[id[1:N]]")) |>
    paste0(collapse = " + ")

  line <- glue::glue("mu[1:N] <-
                    LINPRED(~
                     (x[id[1:N]]) +
                     {expanded_covariates} +
                     (t[1:N]) +
                     (t[1:N]:x[id[1:N]]),
                     coefPrefix=beta_,
                     priors = main_model_priors)")

  return(rlang::parse_expr(line))
}

###############################################################################
create_imputation_model_linpred <- function(imputation_model_covariates,
                                            imputation_model_priors){

  expanded_covariates <-
    purrr::map_chr(imputation_model_covariates,
                   \(x) glue::glue("{x}[1:G]")) |>
    paste0(collapse = " + ")

  line <- glue::glue("eta[1:G] <-
                     LINPRED(~ {expanded_covariates},
                     coefPrefix=gamma_,
                     priors = imputation_model_priors)")
  return(rlang::parse_expr(line))
}

###############################################################################
create_imputation_model_distribution <-
  function(imputation_model_distribution =
             c(
               "normal",
               "binomial",
               "beta_binomial",
               "negative_binomial",
               "poisson",
               "probit_bernoulli"
             )) {
    imputation_model_distribution <- match.arg(imputation_model_distribution)

    switch(imputation_model_distribution,
      normal = {
        line <-
          "
           {
             x[1:G] ~ FORLOOP(dnorm(eta[1:G], sd = sigma_imputation))
             sigma_imputation ~ dexp(rate = 0.1)
           }
           "
      },
      negative_binomial = {
        # https://georgederpa.github.io/teaching/countModels.html
        line <-
          "
           {
             mu_2[1:G] <- FORLOOP(exp(eta[1:G]))
             p[1:G] <- FORLOOP(r / (r + mu_2[1:G]))
             x[1:G] ~ FORLOOP(dnegbin(prob = p[1:G], size = r))
             r ~ dexp(rate = 0.1)
           }
           "
      },
      poisson = {
        line <-
          "
           {
             x[1:G] ~ FORLOOP(dpois(lambda = exp(eta[1:G])))
           }
           "
      },
      binomial = {
        line <-
          "
           {
             x[1:G] ~ FORLOOP(dbinom(prob = expit(eta[1:G]), size = x_size))
           }
           "
      },
      beta_binomial = {
        # https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html#the-beta-binomial-distribution
        line <-
          "
           {

           for(g in 1:G) {
                p[g] <- expit(eta[g])
                s1[g] <- phi * p[g]
                s2[g]  <- phi * (1 - p[g])

                x[g] ~ dBetaBinom_One(N = x_size,
                                      shape1 = s1[g],
                                      shape2 = s2[g])
            }

              phi ~ dexp(rate = 1)

           }
           "
      }
    )

    return(rlang::parse_expr(line))
  }

create_random_effects_code <-
  function(correlated_random_effects = TRUE) {
  if (correlated_random_effects) {
    line <- "
    {
      # Random effects standard deviations
      sigma_u0 ~ dexp(rate = 0.1)  # SD for random intercepts
      sigma_u1 ~ dexp(rate = 0.1)  # SD for random slopes

      # Correlation between random intercept and slope
      rho ~ dunif(-0.99, 0.99)

      # Non-centered random effects (raw parameters)
      for (j in 1:G) {
        u0_raw[j] ~ dnorm(0, sd = 1)  # Raw random intercept
        u1_raw[j] ~ dnorm(0, sd = 1)  # Raw random slope
      }

      # Transform to centered parameterization with correlation
      for (j in 1:G) {
        u0[j] <- sigma_u0 * u0_raw[j]
        u1[j] <- sigma_u1 * (rho * u0_raw[j] + sqrt(1 - rho^2) * u1_raw[j])
      }
    }
    "
  } else {
    line <- "
    {
      # Random effects standard deviations
      sigma_u0 ~ dexp(rate = 0.1)  # SD for random intercepts
      sigma_u1 ~ dexp(rate = 0.1)  # SD for random slopes

      # Non-centered random effects (raw parameters)
      for (j in 1:G) {
        u0_raw[j] ~ dnorm(0, sd = 1)  # Raw random intercept
        u1_raw[j] ~ dnorm(0, sd = 1)  # Raw random slope
      }

      # Transform to centered parameterization without correlation
      for (j in 1:G) {
        u0[j] <- sigma_u0 * u0_raw[j]
        u1[j] <- sigma_u1 * u1_raw[j]
      }
    }
    "
  }

  return(rlang::parse_expr(line))
}

build_nimble_code <- function(main_model_covariates,
                              imputation_model_covariates,
                              imputation_model_distribution,
                              correlated_random_effects = TRUE,
                              main_model_priors,
                              imputation_model_priors) {

  modelCode <- eval(bquote(
    nimbleCode({
      .(create_main_model_linpred(main_model_covariates,
                                  main_model_priors))

      .(create_random_effects_code(correlated_random_effects))

      for (i in 1:N) {
        mu_total[i] <- mu[i] + u0[id[i]] + u1[id[i]] * t[i]
      }

      y[1:N] ~ FORLOOP(dnorm(mu_total[1:N], sd = sigma_residual))
      sigma_residual ~ dexp(rate = 0.1)

      .(create_imputation_model_linpred(imputation_model_covariates,
                                        imputation_model_priors))
      .(create_imputation_model_distribution(imputation_model_distribution))
    })
  ), list(main_model_covariates = main_model_covariates,
          imputation_model_covariates = imputation_model_covariates,
          imputation_model_distribution = imputation_model_distribution,
          correlated_random_effects = correlated_random_effects,
          main_model_priors = main_model_priors,
          imputation_model_priors = imputation_model_priors))

  return(modelCode)
}

#' Fit a Bayesian two-stage model using NIMBLE
#'
#' Fits a mixed effects model with imputation using NIMBLE for MCMC sampling.
#'
#' @param data A data frame containing the outcome and covariates. Must include
#'   columns: y (outcome), t (time), x (exposure), and id (subject identifier).
#' @param main_model_covariates Character vector of covariate names for the main model
#' @param imputation_model_covariates Character vector of covariate names for the imputation model
#' @param imputation_model_distribution Distribution for the imputation model.
#'   One of: "normal", "binomial", "beta_binomial", "poisson", "negative_binomial"
#' @param correlated_random_effects Logical; if TRUE, random intercepts and slopes
#'   are correlated (default: TRUE)
#' @param nchains Number of MCMC chains (default: 4)
#' @param niter Number of MCMC iterations per chain (default: 10000)
#' @param nburnin Number of burn-in iterations to discard (default: 2000)
#' @param x_size Integer vector of trial sizes; required for binomial or beta_binomial
#'   distributions (default: NULL)
#' @param print_summary Logical; if TRUE, prints MCMC summary (default: FALSE)
#' @param print_code Logical; if TRUE, prints the NIMBLE model code (default: FALSE)
#' @return A list containing MCMC samples and optionally WAIC
#' @export
fit_model <- function(data,
                      main_model_covariates,
                      imputation_model_covariates,
                      imputation_model_distribution,
                      correlated_random_effects = TRUE,
                      nchains = 4,
                      niter = 10000,
                      nburnin = 2000,
                      x_size = NULL,
                      print_summary = FALSE,
                      print_code = FALSE
){

  if (imputation_model_distribution %in% c("binomial", "beta_binomial")) {
    stopifnot("Must specify x_size for binomial or beta binomial distributions" = !is.null(x_size))
  }

  check_cols(data, c("y", "t", "x", "id"))

  # Dataframe with one row per ID
  G_df <- dplyr::distinct(data, id, .keep_all = TRUE)

  constants <- list()

  # Scalars
  constants[["N"]] <- nrow(data)
  constants[["G"]] <- nrow(G_df)

  # Length N
  constants[["y"]] <- data[["y"]]
  constants[["t"]] <- data[["t"]]
  constants[["id"]] <- data[["id"]]
  constants[["id_factor"]] <- as.factor(data[["id"]])

  # Size for binomial and beta binomial
  if (!is.null(x_size)) {
    stopifnot("x_size must be integer" = is.integer(x_size))
    constants[["x_size"]] <- as.integer(x_size)
  }

  # Length G
  constants[["x"]] <- G_df[["x"]]

  # All covariate must be time-invariant
  # i.e., baseline
  covariates <- union(main_model_covariates,
                      imputation_model_covariates)

  for (covariate in covariates) {
    constants[[covariate]] <- G_df[[covariate]]
  }

  main_model_priors <- setPriors(
    intercept = quote(dnorm(0, sd = 2)),
    coefficient = quote(dnorm(0, sd = 2))
  )

  if (imputation_model_distribution == "normal") {
    imputation_model_priors <- setPriors(
      intercept = quote(dnorm(0, sd = 2)),
      coefficient = quote(dnorm(0, sd = 2))
    )
  } else if (imputation_model_distribution %in% c("binomial", "beta_binomial")) {
    imputation_model_priors <- setPriors(
      intercept = quote(dnorm(0, sd = 2.5)),
      coefficient = quote(dnorm(0, sd = 2.5))
    )
  } else if (imputation_model_distribution %in% c("poisson", "negative_binomial")) {
    imputation_model_priors <- setPriors(
      intercept = quote(dnorm(0, sd = 2.5)),
      coefficient = quote(dnorm(0, sd = 2.5))
    )
  } else {
    # Default for any other distributions
    imputation_model_priors <- setPriors(
      intercept = quote(dnorm(0, sd = 100)),
      coefficient = quote(dnorm(0, sd = 100))
    )
  }

  code <- build_nimble_code(
    main_model_covariates = main_model_covariates,
    imputation_model_covariates = imputation_model_covariates,
    imputation_model_distribution = imputation_model_distribution,
    correlated_random_effects = correlated_random_effects,
    main_model_priors = main_model_priors,
    imputation_model_priors = imputation_model_priors
  )

  mod <- nimbleModel(code = code,
                     constants = constants,
                     buildDerivs = FALSE)

  # model_code <- mod$getCode()
  # print(model_code)

  vars <- mod$getVarNames()

  # Conditionally include rho in monitoring based on correlation toggle
  if (correlated_random_effects) {
    mon  <- vars[ grepl("^(beta_|gamma_|sigma_)", vars) | vars == "rho" ]
  } else {
    mon  <- vars[ grepl("^(beta_|gamma_|sigma_)", vars) ]
  }

  conf <- configureMCMC(mod,
                        monitors = mon,
                        print = FALSE,
                        enableWAIC = use_waic)

  ##############################################################################

  ##############################################################################

  # conf$printSamplers()

  mcmc  <- buildMCMC(conf)

  Cmod  <- compileNimble(mod)

  Cmcmc <- compileNimble(mcmc,
                         project = Cmod)

  samples <- runMCMC(Cmcmc,
                     nchains = nchains,
                     niter = niter,
                     nburnin = nburnin,
                     WAIC = use_waic,
                     samplesAsCodaMCMC = TRUE)

  if (print_code) {
    cat("CODE:----------------------------------------------------------------")
    print(mod$getCode())
    cat("---------------------------------------------------------------------")
    conf$printSamplers()
    cat("---------------------------------------------------------------------")
  }

  if (print_summary) {
    print(mcmc_summary(samples$samples, dataset_id = "print"))

    if (use_waic) print(samples$WAIC)
  }

  out <- samples$samples

  if (use_waic) {
    attr(out, "WAIC") <- samples$WAIC$WAIC
  }

  return(out)
}
