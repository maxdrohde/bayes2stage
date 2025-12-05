#' Generate correlated random intercept / slope longitudinal data
#'
#' Generates simulated longitudinal data with correlated random intercepts
#' and slopes for testing two-stage design methods.
#'
#' @param N Number of subjects (default: 2000)
#' @param M Number of time points per subject (default: 5)
#' @param alpha_main Intercept for the main model (default: 1)
#' @param beta_x Coefficient for x in the main model (default: 1)
#' @param beta_z Coefficient for z in the main model (default: 1)
#' @param beta_t Coefficient for time in the main model (default: 2)
#' @param beta_t_x_interaction Coefficient for time-x interaction (default: 0.3)
#' @param beta_t_z_interaction Coefficient for time-z interaction (default: 0)
#' @param error_sd Standard deviation of residual error (default: 4)
#' @param x_dist Distribution for x: "normal", "poisson", "binomial",
#'   "negative_binomial", or "beta_binomial" (default: "normal")
#' @param x_size Size parameter for binomial/beta_binomial distributions
#' @param x_disp_param Dispersion parameter for negative_binomial/beta_binomial
#' @param rand_intercept_sd Standard deviation of random intercepts (default: 3)
#' @param rand_slope_sd Standard deviation of random slopes (default: 1)
#' @param rand_eff_corr Correlation between random intercepts and slopes (default: 0)
#' @param gamma0 Intercept for the imputation model (default: 1)
#' @param gamma1 Linear coefficient for z in the imputation model (default: 1)
#' @param gamma2 Quadratic coefficient for z in the imputation model (default: 0)
#' @param gamma_sd Standard deviation for the imputation model (normal distribution only) (default: 2)
#' @return A data frame with simulated longitudinal data
#' @export
generate_data <-
  function(N = 2000,
           M = 5,
           alpha_main = 1, beta_x = 1, beta_z = 1,
           beta_t = 2,
           beta_t_x_interaction = 0.3,
           beta_t_z_interaction = 0,
           error_sd = 4,
           x_dist = c("normal",
                      "poisson",
                      "binomial",
                      "negative_binomial",
                      "beta_binomial"),
           x_size = NULL,
           x_disp_param = NULL,
           rand_intercept_sd = 3, rand_slope_sd = 1, rand_eff_corr = 0,
           gamma0 = 1, gamma1 = 1, gamma2 = 0, gamma_sd = 2) {

    # Checks ###################################################################
    stopifnot("X distribution not supported" = x_dist %in% c("normal",
                                                             "poisson",
                                                             "binomial",
                                                             "negative_binomial",
                                                             "beta_binomial"))

    if (x_dist %in% c("binomial", "beta_binomial")) {
      stopifnot("Must specify x_size for binomial or beta binomial distributions" = !is.null(x_size))
    }

    if (x_dist %in% c("negative_binomial", "beta_binomial")) {
      stopifnot("Must specify x_disp_param for negative binomial and beta binomial distribution" = !is.null(x_disp_param))
    }
    ############################################################################

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
      # Z is standard normal
      z <- stats::rnorm(n = 1, mean = 0, sd = 1)

      eta <-
        gamma0 +
        gamma1*z +
        gamma2*z^2

      x <-
        switch(
          x_dist,
          normal = stats::rnorm(n = 1, mean = eta, sd = gamma_sd),
          poisson = stats::rpois(n = 1, lambda = exp(eta)),
          binomial = stats::rbinom(n = 1, size = x_size, prob = stats::plogis(eta)),
          negative_binomial = stats::rnbinom(n = 1, mu = exp(eta), size = x_disp_param),
          beta_binomial = extraDistr::rbbinom(n = 1,
                                              size = x_size,
                                              alpha = x_disp_param * stats::plogis(eta),
                                              beta = x_disp_param * (1 - stats::plogis(eta)))
        )

      # Generate outcome
      y <-
        alpha_main +
        beta_x * x +
        beta_z * z +
        (beta_t + rand_effs[[2]]) * t +
        (beta_t_x_interaction * t * x) +
        (beta_t_z_interaction * t * z) +
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
