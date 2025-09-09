#' @import nimble
#' @import nimbleEcology
#' @import nimbleMacros
#' @import nimbleHMC

nimbleOptions("enableWAIC" = TRUE)

# Create the NIMBLE model line for the main model
# Allows the user to input any number of covariates

create_main_model_linpred <- function(main_model_covariates){

  expanded_covariates <-
    purrr::map_chr(main_model_covariates, \(x) glue::glue("{x}[id[1:N]]")) |>
    paste0(collapse = " + ")

  #line <- glue::glue("mu[1:N] <- LINPRED(~ (x[id[1:N]]) + {expanded_covariates} + (t[1:N]) + (t[1:N]:x[id[1:N]]) + (t[1:N]|id_factor[1:N]), coefPrefix=beta_, priors = priors)")
  line <- glue::glue("mu[1:N] <- LINPRED(~ x[id[1:N]] + {expanded_covariates} + t[1:N] + t[1:N]:x[id[1:N]] + (t[1:N]|id_factor[1:N]), coefPrefix=beta_, noncentered=TRUE, priors = priors)")

  return(rlang::parse_expr(line))
}

create_imputation_model_linpred <- function(imputation_model_covariates){

  expanded_covariates <-
    purrr::map_chr(imputation_model_covariates, \(x) glue::glue("{x}[1:G]")) |>
    paste0(collapse = " + ")

  line <- glue::glue("eta[1:G] <- LINPRED(~ {expanded_covariates}, coefPrefix=gamma_)")
  return(rlang::parse_expr(line))
}

create_imputation_model_distribution <- function(imputation_model_distribution =
                                                   c("normal",
                                                     "binomial",
                                                     "beta_binomial",
                                                     "negative_binomial",
                                                     "poisson")) {
  imputation_model_distribution <- match.arg(imputation_model_distribution)

  switch(imputation_model_distribution,
         normal = {
           line =
           "
           {
             x[1:G] ~ FORLOOP(dnorm(eta[1:G], sd = sigma_imputation))
             sigma_imputation ~ dexp(rate = 1)
           }
           "
         },
         negative_binomial = {
           # https://georgederpa.github.io/teaching/countModels.html
           line =
             "
           {
             mu_2[1:G] <- FORLOOP(exp(eta[1:G]))
             p[1:G] <- FORLOOP(r / (r + mu_2[1:G]))
             x[1:G] ~ FORLOOP(dnegbin(prob = p[1:G], size = r))
             r ~ dexp(rate = 1)
           }
           "
         },
         poisson = {
           line =
             "
           {
             x[1:G] ~ FORLOOP(dpois(lambda = exp(eta[1:G])))
           }
           "
         },
         binomial = {
           line =
             "
           {
             x[1:G] ~ FORLOOP(dbinom(prob = expit(eta[1:G]), size = x_size))
           }
           "
         },
         beta_binomial = {
        # https://cran.r-project.org/web/packages/brms/vignettes/brms_customfamilies.html#the-beta-binomial-distribution
           line =
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

              phi ~ dexp(rate = lambda)

           }
           "
         }
  )

  return(rlang::parse_expr(line))
}

build_nimble_code <- function(main_model_covariates,
                              imputation_model_covariates,
                              imputation_model_distribution,
                              priors = priors) {

  modelCode <- eval(bquote(
    nimbleCode({
      .(create_main_model_linpred(main_model_covariates))

      y[1:N] ~ FORLOOP(dnorm(mu[1:N], sd = sigma_residual))
      sigma_residual ~ dexp(rate = 1)

      .(create_imputation_model_linpred(imputation_model_covariates))
      .(create_imputation_model_distribution(imputation_model_distribution))
    })
  ), list(main_model_covariates = main_model_covariates,
          imputation_model_covariates = imputation_model_covariates,
          imputation_model_distribution = imputation_model_distribution))

  return(modelCode)
}

#' @export
fit_model <- function(df,
                      main_model_covariates,
                      imputation_model_covariates,
                      imputation_model_distribution,
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

  check_cols(df, c("y", "t", "x", "id"))

  # Dataframe with one row per ID
  G_df <- dplyr::distinct(df, id, .keep_all = TRUE)

  constants <- list()

  # Scalars
  constants[["N"]] <- nrow(df)
  constants[["G"]] <- nrow(G_df)

  # Length N
  constants[["y"]] <- df[["y"]]
  constants[["t"]] <- df[["t"]]
  constants[["id"]] <- df[["id"]]
  constants[["id_factor"]] <- as.factor(df[["id"]])

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

  priors <- setPriors(
    intercept = quote(dnorm(0, sd = 10)),
    coefficient = quote(dnorm(0, sd = 10)),
    sd = quote(dexp(0.5)),
    lkjShape = 1
  )

  code <- build_nimble_code(
    main_model_covariates       = main_model_covariates,
    imputation_model_covariates = imputation_model_covariates,
    imputation_model_distribution = imputation_model_distribution,
    priors = priors
  )

  mod <- nimbleModel(code = code,
                     constants = constants,
                     buildDerivs = TRUE)

  model_code <- mod$getCode()
  print(model_code)

  conf <- configureMCMC(mod,
                       print = FALSE,
                       enableWAIC = TRUE)

  ##############################################################################

  ##############################################################################

  conf$printSamplers()

  mcmc  <- buildMCMC(conf)

  Cmod  <- compileNimble(mod)

  Cmcmc <- compileNimble(mcmc,
                         project = Cmod)

  samples <- runMCMC(Cmcmc,
                     nchains = nchains,
                     niter = niter,
                     nburnin = nburnin,
                     WAIC = TRUE,
                     samplesAsCodaMCMC = TRUE)

  if (print_code) {
    cat("CODE:----------------------------------------------------------------")
    print(mod$getCode())
    cat("---------------------------------------------------------------------")
  }

  if (print_summary) {
    print(mcmc_summary(samples$samples, dataset_id = "print"))
    print(samples$WAIC)
  }

  out <- samples$samples
  attr(out, "WAIC") <- samples$WAIC$WAIC

  return(out)
}
