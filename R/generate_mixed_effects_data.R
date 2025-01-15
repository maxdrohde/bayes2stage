generate_one_subject <- function(M,
                                 id,
                                 alpha_main,
                                 beta_x_e,
                                 beta_x_z,
                                 beta_t,
                                 beta_t_xe_interaction,
                                 error_param,
                                 error_type,
                                 rand_intercept_sd,
                                 rand_slope_sd,
                                 rand_eff_corr,
                                 x_cov,
                                 xe_var,
                                 xz_var){

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

  X <- MASS::mvrnorm(n = 1,
                     mu = c(0,0),
                     Sigma = matrix(c(xe_var,
                                      x_cov,
                                      x_cov,
                                      xz_var),
                                    nrow = 2))

  # x_e is the "expensive" covariate
  x_e <- X[[1]]
  x_z <- X[[2]]

  # Generate outcome
  y <-
    alpha_main +
    beta_x_e * x_e +
    beta_x_z * x_z +
    (beta_t + rand_effs[[2]]) * t +
    (beta_t_xe_interaction * t * x_e) +
    rand_effs[[1]]

  stopifnot("Error type not supported" = error_type %in% c("normal", "uniform", "exponential"))

  y <-
    y +
    switch(
      error_type,
      normal = stats::rnorm(n = M, mean = 0, sd = error_param),
      uniform = stats::runif(n = M, min = -error_param, max = error_param),
      exponential = stats::rexp(n = M, rate = error_param) - (1 / error_param)
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
#' @param error_param For "normal": SD of error term. For "uniform", the range is -error_param to error_param. For "exponential", the rate.
#' @param rand_intercept_sd SD of random intercepts
#' @param rand_slope_sd SD of random slopes
#' @param rand_eff_corr Correlation between random effects
#' @param x_cov Covariance of the covariates
#' @param xe_var Variance of x_e
#' @param xz_var Variance of x_z
#' @return A simulated mixed-effects dataset
#' @export
generate_mixed_effects_data <- function(N = 1000,
                                        Ms = 5,
                                        alpha_main = 1,
                                        beta_x_e = 2,
                                        beta_x_z = 0.5,
                                        beta_t = 2,
                                        beta_t_xe_interaction = 1,
                                        error_param = 1,
                                        error_type = "normal",
                                        rand_intercept_sd = 3,
                                        rand_slope_sd = 3,
                                        rand_eff_corr = 0.4,
                                        x_cov = 0.6,
                                        xe_var = 1,
                                        xz_var = 1){

  ##############################################################################
  # Checks
  stopifnot("N must be divisible by the length of Ms" = (N / length(Ms)) == as.integer(N / length(Ms)))

  ##############################################################################


  df <-
    purrr::map2(1:N,
                rep(Ms, each = N / length(Ms)),
                ~generate_one_subject(M = .y,
                                      id = .x,
                                      alpha_main = alpha_main,
                                      beta_x_e = beta_x_e,
                                      beta_x_z = beta_x_z,
                                      beta_t = beta_t,
                                      beta_t_xe_interaction = beta_t_xe_interaction,
                                      error_type = error_type,
                                      error_param = error_param,
                                      rand_intercept_sd = rand_intercept_sd,
                                      rand_slope_sd = rand_slope_sd,
                                      rand_eff_corr = rand_eff_corr,
                                      x_cov = x_cov,
                                      xe_var = xe_var,
                                      xz_var = xz_var),
                .progress = TRUE) |>
    purrr::list_rbind()

  return(df)
}
