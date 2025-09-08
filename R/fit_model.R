#' @import nimble
#' @import nimbleEcology
#' @import nimbleMacros

nimbleOptions("enableWAIC" = TRUE)

# Create the NIMBLE model line for the main model
# Allows the user to input any number of covariates

create_main_model_linpred <- function(main_model_covariates){

  expanded_covariates <-
    purrr::map_chr(main_model_covariates, \(x) glue::glue("{x}[id[1:N]]")) |>
    paste0(collapse = " + ")

  line <- glue::glue("mu[1:N] <- LINPRED(~ x[id[1:N]] + {expanded_covariates} + t[1:N] + t[1:N]:x[id[1:N]], coefPrefix=beta_)")

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
             sigma_imputation ~ T(dnorm(0, sd = 3), 0, )
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
             r ~ T(dnorm(0, sd = 3), 0, )
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

              phi ~ T(dnorm(0, sd = 2), 0, )

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




      # --- non-centered random effects ---
      sd_id_factor     ~ T(dnorm(0, sd = 3), 0, )
      sd_t_id_factor   ~ T(dnorm(0, sd = 3), 0, )
      re_sds_id_factor[1] <- sd_id_factor
      re_sds_id_factor[2] <- sd_t_id_factor

      Ustar_id_factor[1:2, 1:2] ~ dlkj_corr_cholesky(1, 2)

      U_id_factor[1:2, 1:2] <- uppertri_mult_diag(
        Ustar_id_factor[1:2, 1:2],
        re_sds_id_factor[1:2]
      )

      Ut_id_factor[1:2, 1:2] <- t(U_id_factor[1:2, 1:2])

      for (j in 1:G) {
        z_id_factor[j, 1] ~ dnorm(0, 1)
        z_id_factor[j, 2] ~ dnorm(0, 1)

        beta_id_factor[j]     <- inprod(Ut_id_factor[1:2, 1], z_id_factor[j, 1:2])
        beta_t_id_factor[j]   <- inprod(Ut_id_factor[1:2, 2], z_id_factor[j, 1:2])
      }

      rho_id_factor <- inprod(Ustar_id_factor[1:2, 1],
                              Ustar_id_factor[1:2, 2])
      # --- end non-centered random effects ---

      for (i in 1:N) {
        mu_total[i] <- mu[i] +
          beta_id_factor[id[i]] +
          beta_t_id_factor[id[i]] * t[i]
      }

      y[1:N] ~ FORLOOP(dnorm(mu_total[1:N], sd = sigma_residual))
      sigma_residual ~ T(dnorm(0, sd = 3), 0, )

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

  code <- build_nimble_code(
    main_model_covariates       = main_model_covariates,
    imputation_model_covariates = imputation_model_covariates,
    imputation_model_distribution = imputation_model_distribution
  )

  mod <- nimbleModel(code = code,
                     constants = constants)

  conf <- configureMCMC(mod,
                        print = FALSE,
                        enableWAIC = TRUE)

  ##############################################################################

  conf$removeSamplers(c("sd_id_factor", "sd_t_id_factor"))
  conf$addSampler(
    target  = c("sd_id_factor","sd_t_id_factor"),
    type    = "RW_block",
    control = list(adaptive = TRUE, adaptInterval = 50)
  )

  conf$resetMonitors()
  all_nodes  <- mod$getNodeNames(stochOnly = FALSE, includeData = FALSE)
  selected_nodes <- grep("\\[", all_nodes, value = TRUE, invert = TRUE)
  ustar_nodes <- grep("Ustar|rho", all_nodes, value = TRUE)

  keep <- union(selected_nodes, ustar_nodes)
  conf$addMonitors(keep)

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
