#' Generate correlated random intercept / slope longitudinal data
#' to test two-stage design methods. In the form of Schildcrout (2019)
#' supplement
#' @param N Integer. Total number of subjects in the cohort. Default is 2000.
#' @param Ms Integer vector. Possible numbers of observations per subject,
#'   sampled uniformly. Default is c(4, 5, 6).
#' @param alpha_main Numeric. Fixed intercept (beta_0). Default is 75.
#' @param beta_x Numeric. Main effect of binary exposure X (beta_s). Default is -0.5.
#' @param beta_z Numeric. Effect of continuous confounder Z (beta_c). Default is -2.
#' @param beta_t Numeric. Fixed time effect (beta_t). Default is -1.
#' @param beta_t_x_interaction Numeric. Interaction between time and exposure
#'   (beta_st). Default is -0.5.
#' @param beta_t_z_interaction Numeric. Interaction between time and confounder.
#'   Default is 0.
#' @param error_sd Numeric. Standard deviation of measurement error (sigma_e).
#'   Default is sqrt(12.25) = 3.5.
#' @param x_prevalence Numeric. Prevalence of binary exposure X. Must be between
#'   0 and 1. Default is 0.25.
#' @param rand_intercept_sd Numeric. Standard deviation of random intercept (sigma_0).
#'   Default is sqrt(81) = 9.
#' @param rand_slope_sd Numeric. Standard deviation of random slope (sigma_1).
#'   Default is sqrt(1.56) = 1.25.
#' @param rand_eff_corr Numeric. Correlation between random intercept and slope
#'   (rho). Must be between -1 and 1. Default is 0.
#' @param gamma0 Numeric. Intercept for confounder Z generation. Default is 0.25.
#' @param gamma1 Numeric. Linear coefficient for X in Z generation. Default is 0.5.
#' @param gamma2 Numeric. Quadratic coefficient for X in Z generation. Default is 0.
#' @param gamma_sd Numeric. Standard deviation of Z residual error. Default is 1.
#' @return dataset
#' @export
generate_data_schild_2019_supp <-
  function(N,
           Ms,
           alpha_main,
           beta_x,
           beta_z,
           beta_t,
           beta_t_x_interaction,
           beta_t_z_interaction,
           error_sd,
           x_prevalence,
           rand_intercept_sd,
           rand_slope_sd,
           rand_eff_corr,
           gamma0,
           gamma1,
           gamma2,
           gamma_sd) {

    # List to store the data frame for each subject
    records <- vector(mode = "list", length = N)

    # Loop over each subject and generate a longitudinal record
    for (i in 1:N) {

      # Generate subject specific random effects
      sds <- diag(c(rand_intercept_sd,
                    rand_slope_sd))

      corr <- matrix(c(1, rand_eff_corr,
                       rand_eff_corr, 1),
                     nrow = 2L)

      sigma <- sds %*% corr %*% sds

      rand_effs <-
        MASS::mvrnorm(n = 1L,
                      mu = c(0, 0),
                      Sigma = sigma)

      # Generate covariate values
      M <- sample(Ms, 1L)
      t <- 0L:(M - 1L)

      # X is binary with p = `x_prevalence`
      x <- stats::rbinom(n = 1L,
                         size = 1L,
                         prob = x_prevalence)

      z <- gamma0 +
           gamma1*x +
           gamma2*x^2 +
           stats::rnorm(n = 1L, mean = 0, sd = gamma_sd)

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
    df <- data.table::rbindlist(records) |>
        as.data.frame()

    return(df)
  }
