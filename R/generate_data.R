#' Generate correlated random intercept / slope longitudinal data
#'
#' Generates simulated longitudinal data with correlated random intercepts
#' and slopes for testing two-stage design methods.
#'
#' @param N Number of subjects (default: 2000)
#' @param M Number of time points per subject. Can be a single integer for
#'   fixed observations per subject, or a vector to sample from for variable
#'   observations (default: 5)
#' @param alpha_main Intercept for the main model (default: 1)
#' @param beta_x Coefficient for x in the main model (default: 1)
#' @param beta_z Coefficient for z in the main model (default: 1)
#' @param beta_t Coefficient for time in the main model (default: 2)
#' @param beta_t_x_interaction Coefficient for time-x interaction (default: 0.3)
#' @param beta_t_z_interaction Coefficient for time-z interaction (default: 0)
#' @param error_sd Standard deviation of residual error (default: 4)
#' @param direction Direction of covariate dependency. "x_given_z" generates
#'   z ~ N(0,1) first, then x | z. "z_given_x" generates x ~ Bernoulli first,
#'   then z | x (default: "x_given_z")
#' @param x_dist Distribution for x: "normal", "bernoulli", "poisson",
#'   "negative_binomial", or "beta_binomial". When direction = "z_given_x",
#'   only "bernoulli" is supported (default: "normal")
#' @param x_size Size parameter for binomial/beta_binomial distributions
#' @param x_disp_param Dispersion parameter for negative_binomial/beta_binomial
#' @param x_prevalence Prevalence of binary x when direction = "z_given_x"
#'   (default: 0.25)
#' @param rand_intercept_sd Standard deviation of random intercepts (default: 3)
#' @param rand_slope_sd Standard deviation of random slopes (default: 1)
#' @param rand_eff_corr Correlation between random intercepts and slopes (default: 0)
#' @param gamma0 Intercept for the covariate model (default: 1)
#' @param gamma1 Linear coefficient in the covariate model (default: 1)
#' @param gamma2 Quadratic coefficient in the covariate model (default: 0)
#' @param gamma_sd Standard deviation for normal distribution in covariate
#'   model (default: 2)
#' @param time_scale How to scale the time variable. "normalized" scales time
#'   to the 0-1 interval (default). "integer" uses integer visit numbers
#'   0, 1, 2, ..., M-1 (matches Schildcrout 2019 simulation setup).
#' @return A data frame with simulated longitudinal data
#' @export
generate_data <-
    function(N = 2000L,
             M = 5L,
             alpha_main = 1, beta_x = 1, beta_z = 1,
             beta_t = 2,
             beta_t_x_interaction = 0.3,
             beta_t_z_interaction = 0,
             error_sd = 4,
             direction = c("x_given_z", "z_given_x"),
             x_dist = c("normal",
                        "bernoulli",
                        "binomial",
                        "poisson",
                        "negative_binomial",
                        "beta_binomial"),
             x_size = NULL,
             x_disp_param = NULL,
             x_prevalence = 0.25,
             rand_intercept_sd = 3, rand_slope_sd = 1, rand_eff_corr = 0,
             gamma0 = 1, gamma1 = 1, gamma2 = 0, gamma_sd = 2,
             time_scale = c("normalized", "integer")) {

    # Argument matching
    direction <- match.arg(direction)
    x_dist <- match.arg(x_dist)
    time_scale <- match.arg(time_scale)

    # Validation
    if (direction == "x_given_z") {
        if (x_dist %in% c("binomial", "beta_binomial") && is.null(x_size)) {
            cli::cli_abort("Must specify {.arg x_size} for binomial or beta binomial distributions.")
        }
        if (x_dist %in% c("negative_binomial", "beta_binomial") && is.null(x_disp_param)) {
            cli::cli_abort("Must specify {.arg x_disp_param} for negative binomial and beta binomial distributions.")
        }
    }

    if (direction == "z_given_x" && x_dist != "bernoulli") {
        cli::cli_abort(
            "When {.arg direction} is {.val z_given_x}, {.arg x_dist} must be {.val bernoulli}."
        )
    }

    if (direction == "z_given_x" && (x_prevalence <= 0 || x_prevalence >= 1)) {
        cli::cli_abort("{.arg x_prevalence} must be between 0 and 1 (exclusive).")
    }

    if (any(M < 2L)) {
        cli::cli_abort("{.arg M} values must be >= 2 (at least 2 time points per subject).")
    }

    # List to store the data frame for each subject
    records <- vector(mode = "list", length = N)

    # Pre-compute covariance matrix for random effects (constant)
    sds <- diag(c(rand_intercept_sd, rand_slope_sd))
    corr <- matrix(c(1, rand_eff_corr,
                     rand_eff_corr, 1), nrow = 2L)
    sigma <- sds %*% corr %*% sds

    # Determine if M is fixed or variable
    variable_m <- length(M) > 1L

    # Loop over each subject and generate a longitudinal record
    for (i in seq_len(N)) {

        # Sample number of observations for this subject
        M_i <- if (variable_m) sample(M, 1L) else M

        # Generate subject specific random effects
        rand_effs <-
            MASS::mvrnorm(n = 1L,
                          mu = c(0, 0),
                          Sigma = sigma)

        # Generate time variable
        t <- if (time_scale == "normalized") {
            (0L:(M_i - 1L)) / (M_i - 1L)
        } else {
            0L:(M_i - 1L)
        }

        # Generate covariates based on direction
        if (direction == "x_given_z") {
            # Z first, then X | Z
            z <- stats::rnorm(n = 1L, mean = 0, sd = 1)

            eta <- gamma0 + gamma1 * z + gamma2 * z^2

            x <- switch(
                x_dist,
                normal = stats::rnorm(n = 1L, mean = eta, sd = gamma_sd),
                bernoulli = stats::rbinom(n = 1L, size = 1L, prob = stats::plogis(eta)),
                poisson = stats::rpois(n = 1L, lambda = exp(eta)),
                binomial = stats::rbinom(n = 1L, size = x_size, prob = stats::plogis(eta)),
                negative_binomial = stats::rnbinom(n = 1L, mu = exp(eta), size = x_disp_param),
                beta_binomial = rbbinom(n = 1L,
                                        size = x_size,
                                        alpha = x_disp_param * stats::plogis(eta),
                                        beta = x_disp_param * (1 - stats::plogis(eta)))
            )
        } else {
            # X first (Bernoulli), then Z | X
            x <- stats::rbinom(n = 1L, size = 1L, prob = x_prevalence)

            eta <- gamma0 + gamma1 * x + gamma2 * x^2

            z <- stats::rnorm(n = 1L, mean = eta, sd = gamma_sd)
        }

        # Generate outcome
        y <-
            alpha_main +
            beta_x * x +
            beta_z * z +
            (beta_t + rand_effs[[2]]) * t +
            (beta_t_x_interaction * t * x) +
            (beta_t_z_interaction * t * z) +
            rand_effs[[1]] +
            stats::rnorm(n = M_i, mean = 0, sd = error_sd)

        records[[i]] <-
            data.frame(y = y,
                       t = t,
                       x = x,
                       z = z,
                       rand_int = rand_effs[[1]],
                       rand_slope = rand_effs[[2]],
                       id = i)
    }

    # Merge all the data frames together
    df <- data.table::rbindlist(records) |>
        as.data.frame()

    return(df)
}
