#' Generate correlated random intercept / slope longitudinal data
#' to test two-stage design methods. In the form of Schildcrout (2019)
#' supplement
#' @param N Integer. Total number of subjects in the cohort. Default is 2000.
#' @param Ms Integer vector. Possible numbers of observations per subject,
#'   sampled uniformly. Default is c(4, 5, 6).
#' @param alpha_main Numeric. Fixed intercept (β₀). Default is 75.
#' @param beta_x Numeric. Main effect of binary exposure X (β_s). Default is -0.5.
#' @param beta_z Numeric. Effect of continuous confounder Z (β_c). Default is -2.
#' @param beta_t Numeric. Fixed time effect (β_t). Default is -1.
#' @param beta_t_x_interaction Numeric. Interaction between time and exposure
#'   (β_st). Default is -0.5.
#' @param beta_t_z_interaction Numeric. Interaction between time and confounder.
#'   Default is 0.
#' @param error_sd Numeric. Standard deviation of measurement error (σ_e).
#'   Default is sqrt(12.25) = 3.5.
#' @param x_prevalence Numeric. Prevalence of binary exposure X. Must be between
#'   0 and 1. Default is 0.25.
#' @param rand_intercept_sd Numeric. Standard deviation of random intercept (σ₀).
#'   Default is sqrt(81) = 9.
#' @param rand_slope_sd Numeric. Standard deviation of random slope (σ₁).
#'   Default is sqrt(1.56) = 1.25.
#' @param rand_eff_corr Numeric. Correlation between random intercept and slope
#'   (ρ). Must be between -1 and 1. Default is 0.
#' @param gamma0 Numeric. Intercept for confounder Z generation. Default is 0.25.
#' @param gamma1 Numeric. Linear coefficient for X in Z generation. Default is 0.5.
#' @param gamma2 Numeric. Quadratic coefficient for X in Z generation. Default is 0.
#' @param gamma_sd Numeric. Standard deviation of Z residual error. Default is 1.
#' @return dataset
#' @export
generate_data_schild_2019_supp <-
  function(N = 2000,
           Ms = c(4,5,6),
           alpha_main = 75,
           beta_x = -0.5,
           beta_z = -2,
           beta_t = -1,
           beta_t_x_interaction = -0.5,
           beta_t_z_interaction = 0,
           error_sd = sqrt(12.25),
           x_prevalence = 0.25,
           rand_intercept_sd = sqrt(81), rand_slope_sd = sqrt(1.56), rand_eff_corr = 0,
           gamma0 = 0.25, gamma1 = 0.5, gamma2 = 0, gamma_sd = 1) {

    # List to store the data frame for each subject
    records <- vector(mode = "list", length = N)

    # Loop over each subject and generate a longitudinal record
    for (i in 1:N) {

      # Generate subject specific random effects
      sds <- diag(c(rand_intercept_sd,
                    rand_slope_sd))

      corr <- matrix(c(1, rand_eff_corr,
                       rand_eff_corr, 1),
                     nrow = 2)

      sigma <- sds %*% corr %*% sds

      rand_effs <-
        MASS::mvrnorm(n = 1,
                      mu = c(0,0),
                      Sigma = sigma)

      # Generate covariate values
      M <- sample(Ms, 1)
      t <- 0:(M-1)

      # X is binary with p = `x_prevalence`
      x <- stats::rbinom(n = 1,
                         size = 1,
                         prob = x_prevalence)

      z <- gamma0 +
           gamma1*x +
           gamma2*x^2 +
           stats::rnorm(n = 1, mean = 0, sd = gamma_sd)

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
