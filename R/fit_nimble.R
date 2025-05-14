#' @import nimble
#' @import nimbleEcology
#' @import nimbleMacros

# Create the NIMBLE model line for the main model
# Allows the user to input any number of covariates

create_main_model_linpred <- function(main_model_covariates){

  expanded_covariates <-
    purrr::map_chr(main_model_covariates, \(x) glue::glue("{x}[id[1:N]]")) |>
    paste0(collapse = " + ")

  line <- glue::glue("mu[1:N] <- LINPRED(~ x[id[1:N]] + {expanded_covariates} + t[1:N] + t[1:N]:x[id[1:N]] + (t[1:N]|id_factor[1:N]), coefPrefix=beta_)")
  # line <- glue::glue("mu[1:N] <- LINPRED(~ x[id[1:N]] + {expanded_covariates} + t[1:N] + t[1:N]:x[id[1:N]] + (t[1:N]||id_factor[1:N]), coefPrefix=beta_, noncentered=TRUE)")
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
                                                     "negative_binomial",
                                                     "poisson")) {
  imputation_model_distribution <- match.arg(imputation_model_distribution)

  ## FIXME: ADD THE OTHER DISTRIBUTIONS OTHER THAN NORMAL
  switch(imputation_model_distribution,
         normal = {
           line =
           "
           {
             x[1:G] ~ FORLOOP(dnorm(eta[1:G], sd = sigma_imputation))
             sigma_imputation ~ dunif(0, 100)
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
             r ~ dunif(0, 100)
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
         # FIXME: NEED TO FIX THIS
         beta_binomial = {
           line =
             "
           {
             x[1:G] ~ FORLOOP(dbinom(prob = expit(eta[1:G]), size = x_size))
           }
           "
         }
  )

  return(rlang::parse_expr(line))
}

build_nimble_code <- function(main_model_covariates,
                              imputation_model_covariates,
                              imputation_model_distribution) {

  modelCode <- eval(bquote(
    nimbleCode({
      .(create_main_model_linpred(main_model_covariates))
      y[1:N] ~ FORLOOP(dnorm(mu[1:N], sd = sigma_residual))
      sigma_residual ~ dunif(0, 100)

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
                      nchain = 4,
                      niter = 10000,
                      nburnin = 2000,
                      x_size = NULL
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

  code <- build_nimble_code(
    main_model_covariates       = main_model_covariates,
    imputation_model_covariates = imputation_model_covariates,
    imputation_model_distribution = imputation_model_distribution
  )

  mod <- nimbleModel(code = code, constants = constants)

  samples <- nimbleMCMC(mod,
                        nchain = nchain,
                        niter = niter,
                        nburnin = nburnin,
                        samplesAsCodaMCMC = TRUE)

  return(samples)
}
