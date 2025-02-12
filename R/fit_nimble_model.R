#' Fit the model with NIMBLE
#' @export
fit_nimble_model <-
  function(df,
           niter = 20000,
           nburnin = 5000,
           nchains = 4,
           print_model_details = FALSE) {

    library(nimble)

    ### Define NIMBLE Model
    code <- nimble::nimbleCode({

      # Main model likelihood
      for(i in 1:N) {
        y[i] ~ dnorm(mu[i], sd = sd_y)

        mu[i] <-
          intercept_main +
          betaX * X[id[i]] +
          betaZ * Zrep[i] +
          betaXt * X[id[i]] * t[i] +
          b0[id[i]] +
          (beta_t + b1[id[i]]) * t[i]
      }

      # Imputation model
      for(g in 1:G) {
        X[g] ~ dnorm(eta[g], sd = sd_gamma)
        eta[g] <- gamma0 + gamma1*Z[g] + gamma2*Z[g]^2
      }

      # Random effects model centered parameterization
      # for(g in 1:G) {
      #   b0[g] ~ dnorm(0, sd = sd_b0)  # random intercept
      #   b1[g] ~ dnorm(0, sd = sd_b1)  # random slope
      # }

      # Non-centered parameterization for random effects:
      for(g in 1:G) {
        b0_raw[g] ~ dnorm(0, 1)          # latent standard normal for random intercept
        b0[g] <- b0_raw[g] * sd_b0        # scale to obtain b0

        b1_raw[g] ~ dnorm(0, 1)          # latent standard normal for random slope
        b1[g] <- b1_raw[g] * sd_b1        # scale to obtain b1
      }

      ### Priors
      intercept_main ~ dnorm(0, sd = 100)
      betaX ~ dnorm(0, sd = 10)
      betaZ ~ dnorm(0, sd = 10)
      beta_t ~ dnorm(0, sd = 10)
      betaXt ~ dnorm(0, sd = 10)

      gamma0 ~ dnorm(0, sd = 100)
      gamma1 ~ dnorm(0, sd = 10)
      gamma2 ~ dnorm(0, sd = 10)


      sd_y  ~ dexp(0.1)
      sd_gamma ~ dexp(0.1)

      sd_b0 ~ dexp(0.1)
      sd_b1 ~ dexp(0.1)
    })


    l <- bayes2stage:::format_data_mcmc(df)

    constants <- l[c("N", "G", "id")]

    data <- l[c("y", "t", "X", "Z", "Zrep")]

    inits <- list(sd_y = 1,
                  sd_b0 = 1,
                  sd_b1 = 1,
                  sd_gamma = 1,
                  intercept_main = 0,
                  betaX = 0,
                  betaZ = 0,
                  beta_t = 0,
                  betaXt = 0,
                  gamma0 = 0,
                  gamma1 = 0,
                  gamma2 = 0)

    model <-
      nimble::nimbleModel(code = code,
                          name = "model",
                          constants = constants,
                          data = data,
                          inits = inits)

    monitors <- c("sd_y", "sd_b0", "sd_b1","sd_gamma",
                  "intercept_main", "betaX", "betaZ", "beta_t",
                  "betaXt", "gamma0", "gamma1", "gamma2")

    # Compile model and MCMC
    Cmodel <- compileNimble(model)
    config_Cmodel <- configureMCMC(Cmodel, print = print_model_details)
    config_Cmodel$resetMonitors()
    config_Cmodel$addMonitors(monitors)
    mcmc <- buildMCMC(config_Cmodel)

    # for (param in c("sd_y", "sd_gamma", "sd_b0", "sd_b1")) {
    #   mcmc$removeSampler(param)
    #   mcmcConf$addSampler(target = param, type = "RW", control = list(log = TRUE))
    # }

    Cmcmc <- compileNimble(mcmc)

    # Run MCMC
    samples <-
      runMCMC(mcmc = Cmcmc,
              niter = niter,
              samplesAsCodaMCMC = TRUE,
              nburnin = nburnin,
              nchains = nchains)

    return(samples)
  }
