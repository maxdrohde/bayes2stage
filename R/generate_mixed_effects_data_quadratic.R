generate_one_subject_quadratic <- function(M,
                                 id,
                                 alpha_main,
                                 beta_x_e,
                                 beta_x_z,
                                 beta_t,
                                 beta_t_xe_interaction,
                                 error_sd,
                                 error_type,
                                 rand_intercept_sd,
                                 rand_slope_sd,
                                 rand_eff_corr,
                                 gamma0,
                                 gamma1,
                                 gamma2,
                                 gamma_sd){

  # Generate subject specific random effects
  sds <- diag(c(rand_intercept_sd, rand_slope_sd))
  corr <- matrix(c(1, rand_eff_corr, rand_eff_corr, 1), nrow = 2)
  sigma <- sds %*% corr %*% sds
  rand_effs <-
    MASS::mvrnorm(n = 1,
                  mu = c(0,0),
                  Sigma = sigma)

  # Generate covariate values
  t <- 0:(M-1)

  x_z <- rnorm(n = 1, mean = 0, sd = 1)
  # Quadratic relationship
  x_e <- gamma0 + gamma1*x_z + gamma2*x_z^2 + rnorm(n = 1, mean = 0, sd = gamma_sd)

  # Generate outcome
  y <-
    alpha_main +
    beta_x_e * x_e +
    beta_x_z * x_z +
    (beta_t + rand_effs[[2]]) * t +
    (beta_t_xe_interaction * t * x_e) +
    rand_effs[[1]]

  stopifnot("Error type not supported" = error_type %in% c("normal", "uniform", "exponential"))

  find_unif_min_max <- function(sd) {
    M <- sd * sqrt(3)
    return(c(min = -M, max = M))
  }

  y <-
    y +
    switch(
      error_type,
      normal = stats::rnorm(n = M, mean = 0, sd = error_sd),
      uniform = stats::runif(n = M, min = -error_sd * sqrt(3), max = error_sd * sqrt(3)),
      exponential = stats::rexp(n = M, rate = 1 / error_sd) - (error_sd)
    )

  return(data.frame(y=y,
                    t=t,
                    x_e=x_e,
                    x_z=x_z,
                    rand_int = rand_effs[[1]],
                    rand_slope = rand_effs[[2]],
                    id=id))
}

#' Generate continuous mixed-effects data for the purpose of testing the Bayesian
#' Two-Stage design methods. We generate longitudinal data with random intercepts
#' and slopes. There are two continuous covariate, the expensive covariate
#' x_e and the inexpensive covariate x_z.
#'
#' @param N Number of subjects to generate
#' @param Ms A vector containing the number of timepoints to generate
#' @param alpha_main Intercept
#' @param beta_x_e x_e effect
#' @param beta_x_z x_z effect
#' @param beta_t time effect
#' @param beta_t_xe_interaction time by x_e interaction
#' @param error_sd SD of error term
#' @param rand_intercept_sd SD of random intercepts
#' @param rand_slope_sd SD of random slopes
#' @param rand_eff_corr Correlation between random effects
#' @return A simulated mixed-effects dataset
#' @export
generate_mixed_effects_data_quadratic <- function(N = 1000,
                                        Ms = 5,
                                        alpha_main = 1,
                                        beta_x_e = 2,
                                        beta_x_z = 0.5,
                                        beta_t = 2,
                                        beta_t_xe_interaction = 1,
                                        error_sd = 1,
                                        error_type = "normal",
                                        rand_intercept_sd = 3,
                                        rand_slope_sd = 3,
                                        rand_eff_corr = 0.4,
                                        gamma0 = 1,
                                        gamma1 = 0,
                                        gamma2 = 2,
                                        gamma_sd = 1){

  ##############################################################################
  # Checks
  stopifnot("N must be divisible by the length of Ms" = (N / length(Ms)) == as.integer(N / length(Ms)))

  ##############################################################################


  df <-
    purrr::map2(1:N,
                rep(Ms, each = N / length(Ms)),
                ~generate_one_subject_quadratic(M = .y,
                                      id = .x,
                                      alpha_main = alpha_main,
                                      beta_x_e = beta_x_e,
                                      beta_x_z = beta_x_z,
                                      beta_t = beta_t,
                                      beta_t_xe_interaction = beta_t_xe_interaction,
                                      error_sd = error_sd,
                                      error_type = error_type,
                                      rand_intercept_sd = rand_intercept_sd,
                                      rand_slope_sd = rand_slope_sd,
                                      rand_eff_corr = rand_eff_corr,
                                      gamma0 = gamma0,
                                      gamma1 = gamma1,
                                      gamma2 = gamma2,
                                      gamma_sd = gamma_sd),
                .progress = TRUE) |>
    purrr::list_rbind()

  return(df)
}
