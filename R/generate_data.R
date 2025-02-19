#' Generate correlated random intercept / slope longitudinal data
#' to test two-stage design methods.
#'
#' @param N TO DO
#' @return TO DO
#' @export
generate_data <-
  function(N = 2000,
           M = 5,
           alpha_main = 1, beta_x = 1, beta_z = 1,
           beta_t = 2, beta_t_xe_interaction = 0.3,
           error_sd = 4, error_type = "normal",
           rand_intercept_sd = 3, rand_slope_sd = 1, rand_eff_corr = 0,
           gamma0 = 1, gamma1 = 1, gamma2 = 0, gamma_sd = 2) {

    # Checks
    # stopifnot("N must be divisible by the length of Ms" = (N / length(Ms)) == as.integer(N / length(Ms)))
    stopifnot("Error type not supported" = error_type %in% c("normal", "uniform", "exponential"))

    # List to store the data frame for each subject
    records <- vector(mode = "list", length = N)

    # Loop over each subject and generate a longitudinal record
    for (i in 1:N) {

      # Generate subject specific random effects
      sds <- diag(c(rand_intercept_sd, rand_slope_sd))
      corr <- matrix(c(1, rand_eff_corr,
                       rand_eff_corr, 1), nrow = 2)
      sigma <- sds %*% corr %*% sds
      rand_effs <-
        MASS::mvrnorm(n = 1,
                      mu = c(0,0),
                      Sigma = sigma)

      # Generate covariate values
      t <- (0:(M-1)) / (M-1)
      z <- rnorm(n = 1, mean = 0, sd = 1)
      x <-
        gamma0 +
        gamma1*z +
        gamma2*z^2

      x <-
        x +
        switch(
          error_type,
          normal = stats::rnorm(n = 1, mean = 0, sd = gamma_sd),
          uniform = stats::runif(n = 1, min = -gamma_sd * sqrt(3), max = gamma_sd * sqrt(3)),
          exponential = stats::rexp(n = 1, rate = 1 / gamma_sd) - (gamma_sd)
        )

      # Generate outcome
      y <-
        alpha_main +
        beta_x * x +
        beta_z * z +
        (beta_t + rand_effs[[2]]) * t +
        (beta_t_xe_interaction * t * x) +
        rand_effs[[1]] +
        stats::rnorm(n = M, mean = 0, sd = error_sd)

      records[[i]] <-
        data.frame(y=y,
                   t=t,
                   x=x,
                   z=z,
                   rand_int = rand_effs[[1]],
                   rand_slope = rand_effs[[2]],
                   id=i)
    }

    # Merge all the data frames together
    df <- do.call(rbind, records)

    return(df)
  }
