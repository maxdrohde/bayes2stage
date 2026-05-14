#' @title Posterior Predictive Checks for Imputation Model
#' @name ppc_imputation
#' @description Functions to diagnose imputation model adequacy
NULL

# =============================================================================
# Helper Functions
# =============================================================================

#' Extract observed x and design matrix Z safely
#' @param data Full data frame
#' @param imputation_covariates Character vector of covariate names
#' @return List with x (vector) and Z (matrix)
#' @keywords internal
get_ppc_data <- function(data, imputation_covariates) {
  subject_data <- data |>
    dplyr::summarize(
      x = dplyr::first(x),
      dplyr::across(dplyr::all_of(imputation_covariates), dplyr::first),
      .by = id
    ) |>
    dplyr::filter(!is.na(x)) |>
    dplyr::arrange(id)

  Z <- as.matrix(subject_data[, imputation_covariates, drop = FALSE])

  return(list(
    x = subject_data$x,
    Z = Z,
    id = subject_data$id
  ))
}

#' Extract gamma matrix from posterior draws
#' @param draws Posterior draws data frame
#' @param draw_idx Indices of draws to use
#' @return Matrix with gamma values (n_draws x n_covariates)
#' @keywords internal
extract_gamma_matrix <- function(draws, draw_idx = NULL) {
  if (is.null(draw_idx)) {
    draw_idx <- seq_len(nrow(draws))
  }


  gamma_cols <- grep("^gamma\\[", names(draws), value = TRUE)

  if (length(gamma_cols) == 0) {
    return(matrix(0, nrow = length(draw_idx), ncol = 0))
  }

  gamma_cols <- gamma_cols[order(as.numeric(gsub("gamma\\[|\\]", "", gamma_cols)))]
  gamma_mat <- as.matrix(draws[draw_idx, gamma_cols, drop = FALSE])

  return(gamma_mat)
}

#' Extract posterior means for gamma
#' @param summary_df Summary data frame from fit$summary()
#' @return Vector of gamma posterior means
#' @keywords internal
get_gamma_posterior_mean <- function(summary_df) {
  gamma_rows <- summary_df |>
    dplyr::filter(grepl("^gamma\\[", variable)) |>
    dplyr::arrange(as.numeric(gsub("gamma\\[|\\]", "", variable)))

  if (nrow(gamma_rows) == 0) {
    return(0)
  }

  return(gamma_rows$mean)
}

# =============================================================================
# Core PPC Functions
# =============================================================================

#' Posterior predictive density check for imputation model
#'
#' @param fit CmdStanMCMC fit object from fit_stan_model()
#' @param data Data frame used for fitting
#' @param imputation_covariates Character vector of imputation model covariates
#' @param n_draws Number of posterior draws to visualize
#' @return ggplot object
#' @export
ppc_imputation_density <- function(fit, data, imputation_covariates, n_draws = 100) {
  ppc_data <- get_ppc_data(data, imputation_covariates)
  x_obs <- ppc_data$x
  Z_obs <- ppc_data$Z
  n_obs <- length(x_obs)

  draws <- posterior::as_draws_df(fit$draws())
  draw_idx <- sample(nrow(draws), min(n_draws, nrow(draws)))

  alpha_imp <- draws$alpha_imputation[draw_idx]
  sigma_imp <- draws$sigma_imputation[draw_idx]
  gamma_mat <- extract_gamma_matrix(draws, draw_idx)

  x_rep_list <- lapply(seq_along(draw_idx), \(m) {
    eta <- alpha_imp[m]
    if (ncol(Z_obs) > 0) {
      eta <- eta + as.vector(Z_obs %*% gamma_mat[m, ])
    }
    x_rep <- stats::rnorm(n_obs, mean = eta, sd = sigma_imp[m])
    data.frame(x = x_rep, draw = m, type = "predicted")
  })

  x_rep_df <- dplyr::bind_rows(x_rep_list)
  x_obs_df <- data.frame(x = x_obs, draw = 0, type = "observed")

  p <- ggplot2::ggplot() +
    ggplot2::geom_density(
      data = x_rep_df,
      ggplot2::aes(x = x, group = draw),
      color = "steelblue", alpha = 0.15, linewidth = 0.3
    ) +
    ggplot2::geom_density(
      data = x_obs_df,
      ggplot2::aes(x = x),
      color = "darkred", linewidth = 1.5
    ) +
    ggplot2::labs(
      title = "Posterior Predictive Check: Imputation Model",
      subtitle = "Blue = predicted x|z draws; Red = observed x",
      x = "x", y = "Density"
    ) +
    ggplot2::theme_minimal()

  return(p)
}

#' Posterior predictive check for a summary statistic
#'
#' @param fit CmdStanMCMC fit object
#' @param data Data frame used for fitting
#' @param imputation_covariates Character vector of imputation model covariates
#' @param stat Function to compute statistic
#' @param stat_name Name for display
#' @param n_draws Number of posterior draws
#' @return List with observed value, predicted distribution, p-value, and plot
#' @export
ppc_imputation_stat <- function(
  fit,
  data,
  imputation_covariates,
  stat = mean,
  stat_name = "Mean",
  n_draws = 1000
) {
  ppc_data <- get_ppc_data(data, imputation_covariates)
  x_obs <- ppc_data$x
  Z_obs <- ppc_data$Z
  n_obs <- length(x_obs)

  T_obs <- stat(x_obs)

  draws <- posterior::as_draws_df(fit$draws())
  draw_idx <- sample(nrow(draws), min(n_draws, nrow(draws)))

  alpha_imp <- draws$alpha_imputation[draw_idx]
  sigma_imp <- draws$sigma_imputation[draw_idx]
  gamma_mat <- extract_gamma_matrix(draws, draw_idx)

  T_rep <- sapply(seq_along(draw_idx), \(m) {
    eta <- alpha_imp[m]
    if (ncol(Z_obs) > 0) {
      eta <- eta + as.vector(Z_obs %*% gamma_mat[m, ])
    }
    x_rep <- stats::rnorm(n_obs, mean = eta, sd = sigma_imp[m])
    stat(x_rep)
  })

  p_lower <- mean(T_rep <= T_obs)
  p_upper <- mean(T_rep >= T_obs)
  p_value <- 2 * min(p_lower, p_upper)

  p <- ggplot2::ggplot(data.frame(T = T_rep), ggplot2::aes(x = T)) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = T_obs, color = "darkred", linewidth = 1.5) +
    ggplot2::labs(
      title = sprintf("PPC: %s of x", stat_name),
      subtitle = sprintf("Observed = %.3f, p = %.3f", T_obs, p_value),
      x = stat_name, y = "Count"
    ) +
    ggplot2::theme_minimal()

  return(list(
    observed = T_obs,
    predicted = T_rep,
    p_value = p_value,
    p_lower = p_lower,
    p_upper = p_upper,
    plot = p
  ))
}

#' Check prediction interval coverage
#'
#' @param fit CmdStanMCMC fit object
#' @param data Data frame
#' @param imputation_covariates Character vector of imputation model covariates
#' @param prob Coverage probability
#' @param n_draws Number of posterior draws
#' @return List with coverage proportion and plot
#' @export
ppc_imputation_intervals <- function(
  fit,
  data,
  imputation_covariates,
  prob = 0.95,
  n_draws = 1000
) {
  ppc_data <- get_ppc_data(data, imputation_covariates)
  x_obs <- ppc_data$x
  Z_obs <- ppc_data$Z
  n_obs <- length(x_obs)

  draws <- posterior::as_draws_df(fit$draws())
  draw_idx <- sample(nrow(draws), min(n_draws, nrow(draws)))

  alpha_imp <- draws$alpha_imputation[draw_idx]
  sigma_imp <- draws$sigma_imputation[draw_idx]
  gamma_mat <- extract_gamma_matrix(draws, draw_idx)

  x_rep_mat <- matrix(NA, nrow = length(draw_idx), ncol = n_obs)
  for (m in seq_along(draw_idx)) {
    eta <- alpha_imp[m]
    if (ncol(Z_obs) > 0) {
      eta <- eta + as.vector(Z_obs %*% gamma_mat[m, ])
    }
    x_rep_mat[m, ] <- stats::rnorm(n_obs, mean = eta, sd = sigma_imp[m])
  }

  alpha_level <- (1 - prob) / 2
  lower <- apply(x_rep_mat, 2, stats::quantile, probs = alpha_level)
  upper <- apply(x_rep_mat, 2, stats::quantile, probs = 1 - alpha_level)
  pred_mean <- colMeans(x_rep_mat)

  covered <- (x_obs >= lower) & (x_obs <= upper)
  coverage <- mean(covered)

  sort_order <- order(pred_mean)
  plot_df <- data.frame(
    idx = seq_along(pred_mean),
    x_obs = x_obs[sort_order],
    lower = lower[sort_order],
    upper = upper[sort_order],
    pred_mean = pred_mean[sort_order],
    covered = covered[sort_order]
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = idx)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      fill = "steelblue", alpha = 0.3
    ) +
    ggplot2::geom_line(ggplot2::aes(y = pred_mean), color = "steelblue") +
    ggplot2::geom_point(
      ggplot2::aes(y = x_obs, color = covered),
      size = 2
    ) +
    ggplot2::scale_color_manual(
      values = c("TRUE" = "darkgreen", "FALSE" = "red"),
      labels = c("TRUE" = "Yes", "FALSE" = "No")
    ) +
    ggplot2::labs(
      title = sprintf("Prediction Interval Coverage (%.0f%%)", 100 * prob),
      subtitle = sprintf("Actual coverage: %.1f%%", 100 * coverage),
      x = "Subject (sorted by predicted mean)",
      y = "x",
      color = "Covered"
    ) +
    ggplot2::theme_minimal()

  return(list(
    coverage = coverage,
    expected_coverage = prob,
    n_covered = sum(covered),
    n_total = n_obs,
    plot = p,
    data = plot_df
  ))
}

#' Residual diagnostics for normal imputation model
#'
#' @param fit CmdStanMCMC fit object
#' @param data Data frame
#' @param imputation_covariates Character vector of imputation model covariates
#' @return List with diagnostic plots
#' @export
ppc_imputation_residuals <- function(fit, data, imputation_covariates) {
  ppc_data <- get_ppc_data(data, imputation_covariates)
  x_obs <- ppc_data$x
  Z_obs <- ppc_data$Z

  summary_df <- fit$summary()

  alpha_mean <- summary_df |>
    dplyr::filter(variable == "alpha_imputation") |>
    dplyr::pull(mean)

  sigma_mean <- summary_df |>
    dplyr::filter(variable == "sigma_imputation") |>
    dplyr::pull(mean)

  gamma_mean <- get_gamma_posterior_mean(summary_df)

  fitted_vals <- alpha_mean
  if (ncol(Z_obs) > 0 && length(gamma_mean) == ncol(Z_obs)) {
    fitted_vals <- fitted_vals + as.vector(Z_obs %*% gamma_mean)
  }
  residuals <- x_obs - fitted_vals
  std_residuals <- residuals / sigma_mean

  resid_df <- data.frame(
    x_obs = x_obs,
    fitted = fitted_vals,
    residual = residuals,
    std_residual = std_residuals
  )

  plots <- list()

  plots$vs_fitted <- ggplot2::ggplot(
    resid_df,
    ggplot2::aes(x = fitted, y = residual)
  ) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "blue", formula = y ~ x) +
    ggplot2::labs(
      title = "Residuals vs Fitted",
      x = "Fitted E[x|z]",
      y = "Residual"
    ) +
    ggplot2::theme_minimal()

  # Residuals vs all covariates (faceted)
  if (ncol(Z_obs) > 0 && length(imputation_covariates) > 0) {
    resid_vs_z_long <- purrr::imap_dfr(imputation_covariates, \(cov, i) {
      data.frame(
        covariate = cov,
        z_value = Z_obs[, i],
        residual = residuals
      )
    })

    plots$vs_z <- ggplot2::ggplot(
      resid_vs_z_long,
      ggplot2::aes(x = z_value, y = residual)
    ) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
      ggplot2::geom_smooth(method = "loess", se = TRUE, color = "blue", formula = y ~ x) +
      ggplot2::facet_wrap(~covariate, scales = "free_x") +
      ggplot2::labs(
        title = "Residuals vs Covariates",
        subtitle = "Curvature suggests non-linear z-x relationship",
        x = "Covariate value",
        y = "Residual"
      ) +
      ggplot2::theme_minimal()
  } else {
    plots$vs_z <- NULL
  }

  plots$qq <- ggplot2::ggplot(
    resid_df,
    ggplot2::aes(sample = std_residual)
  ) +
    ggplot2::geom_qq() +
    ggplot2::geom_qq_line(color = "red") +
    ggplot2::labs(
      title = "Q-Q Plot of Standardized Residuals",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    ggplot2::theme_minimal()

  plots$scale_location <- ggplot2::ggplot(
    resid_df,
    ggplot2::aes(x = fitted, y = sqrt(abs(std_residual)))
  ) +
    ggplot2::geom_point(alpha = 0.6) +
    ggplot2::geom_smooth(method = "loess", se = FALSE, color = "blue", formula = y ~ x) +
    ggplot2::labs(
      title = "Scale-Location Plot",
      subtitle = "Trend suggests heteroscedasticity",
      x = "Fitted",
      y = expression(sqrt("|Standardized Residual|"))
    ) +
    ggplot2::theme_minimal()

  plots$histogram <- ggplot2::ggplot(
    resid_df,
    ggplot2::aes(x = std_residual)
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      bins = 30, fill = "steelblue", color = "white"
    ) +
    ggplot2::stat_function(fun = stats::dnorm, color = "red", linewidth = 1) +
    ggplot2::labs(
      title = "Distribution of Standardized Residuals",
      subtitle = "Red = standard normal",
      x = "Standardized Residual",
      y = "Density"
    ) +
    ggplot2::theme_minimal()

  n_test <- min(length(std_residuals), 5000)
  shapiro_result <- stats::shapiro.test(std_residuals[seq_len(n_test)])

  diagnostics <- list(
    mean_residual = mean(residuals),
    sd_residual = stats::sd(residuals),
    shapiro_p = shapiro_result$p.value
  )

  return(list(
    plots = plots,
    diagnostics = diagnostics,
    data = resid_df
  ))
}

# =============================================================================
# Mean Model Diagnostics
# =============================================================================

#' Mean model diagnostics for normal imputation
#'
#' Comprehensive diagnostics to detect misspecification of the conditional mean
#' E\[x|z\] in the normal imputation model. Produces four diagnostic plots that
#' help identify nonlinearity, systematic patterns, and model fit issues.
#'
#' @param fit CmdStanMCMC fit object from fit_stan_model()
#' @param data Data frame used for fitting
#' @param imputation_covariates Character vector of imputation model covariate names
#' @param n_bins Number of bins for binned residual plot (default: 10)
#' @return A list containing:
#'   \item{plots}{List of ggplot objects: resid_vs_z, partial_residuals, obs_vs_pred, binned_residuals}
#'   \item{r_squared}{R-squared of observed vs predicted}
#'   \item{fitted}{Vector of fitted values E\[x|z\]}
#'   \item{residuals}{Vector of residuals (observed - fitted)}
#' @export
ppc_imputation_mean_diagnostics <- function(
    fit,
    data,
    imputation_covariates,
    n_bins = 10
) {
    ppc_data <- get_ppc_data(data, imputation_covariates)
    x_obs <- ppc_data$x
    Z_obs <- ppc_data$Z

    if (length(imputation_covariates) == 0 || ncol(Z_obs) == 0) {
        cli::cli_abort("No imputation covariates provided")
    }

    summary_df <- fit$summary()
    alpha_mean <- summary_df |>
        dplyr::filter(variable == "alpha_imputation") |>
        dplyr::pull(mean)
    gamma_mean <- get_gamma_posterior_mean(summary_df)

    if (length(gamma_mean) == 0) {
        gamma_mean <- rep(0, ncol(Z_obs))
    }

    fitted_vals <- as.vector(alpha_mean + Z_obs %*% gamma_mean)
    residuals <- x_obs - fitted_vals

    plots <- list()

    # 1. Residuals vs all covariates (faceted)
    resid_long <- purrr::imap_dfr(imputation_covariates, \(cov, i) {
        data.frame(
            covariate = cov,
            z_value = Z_obs[, i],
            residual = residuals
        )
    })

    plots$resid_vs_z <- ggplot2::ggplot(
        resid_long,
        ggplot2::aes(x = z_value, y = residual)
    ) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        ggplot2::geom_smooth(method = "loess", se = TRUE, color = "blue",
                             formula = y ~ x) +
        ggplot2::facet_wrap(~covariate, scales = "free_x") +
        ggplot2::labs(
            title = "Residuals vs Covariates",
            subtitle = "Curvature in loess suggests nonlinear relationship",
            x = "Covariate value", y = "Residual"
        ) +
        ggplot2::theme_minimal()

    # 2. Partial residual plots (faceted)
    partial_long <- purrr::imap_dfr(imputation_covariates, \(cov, i) {
        partial_resid <- gamma_mean[i] * Z_obs[, i] + residuals
        data.frame(
            covariate = cov,
            z_value = Z_obs[, i],
            partial_residual = partial_resid
        )
    })

    plots$partial_residuals <- ggplot2::ggplot(
        partial_long,
        ggplot2::aes(x = z_value, y = partial_residual)
    ) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red",
                             formula = y ~ x) +
        ggplot2::geom_smooth(method = "loess", se = TRUE, color = "blue",
                             formula = y ~ x) +
        ggplot2::facet_wrap(~covariate, scales = "free") +
        ggplot2::labs(
            title = "Partial Residual Plots",
            subtitle = "Red = linear fit; Blue = loess (should overlap if linear is correct)",
            x = "Covariate value", y = "Partial residual"
        ) +
        ggplot2::theme_minimal()

    # 3. Observed vs predicted
    r_squared <- stats::cor(x_obs, fitted_vals)^2
    obs_pred_df <- data.frame(observed = x_obs, predicted = fitted_vals)

    plots$obs_vs_pred <- ggplot2::ggplot(
        obs_pred_df,
        ggplot2::aes(x = predicted, y = observed)
    ) +
        ggplot2::geom_point(alpha = 0.5) +
        ggplot2::geom_abline(slope = 1, intercept = 0, color = "red",
                             linetype = "dashed") +
        ggplot2::labs(
            title = "Observed vs Predicted x",
            subtitle = sprintf("R-squared = %.3f", r_squared),
            x = "Predicted E[x|z]", y = "Observed x"
        ) +
        ggplot2::theme_minimal()

    # 4. Binned residuals
    bin_df <- data.frame(fitted = fitted_vals, residual = residuals) |>
        dplyr::mutate(bin = cut(fitted, breaks = n_bins, include.lowest = TRUE)) |>
        dplyr::summarize(
            mean_fitted = mean(fitted),
            mean_resid = mean(residual),
            se_resid = stats::sd(residual) / sqrt(dplyr::n()),
            n = dplyr::n(),
            .by = bin
        )

    plots$binned_residuals <- ggplot2::ggplot(
        bin_df,
        ggplot2::aes(x = mean_fitted, y = mean_resid)
    ) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = -1.96 * se_resid, ymax = 1.96 * se_resid),
            alpha = 0.2, fill = "steelblue"
        ) +
        ggplot2::geom_point(ggplot2::aes(size = n), color = "steelblue") +
        ggplot2::geom_line(color = "steelblue") +
        ggplot2::labs(
            title = "Binned Residuals",
            subtitle = "Points outside bands suggest systematic misfit",
            x = "Mean predicted value (binned)", y = "Mean residual",
            size = "N per bin"
        ) +
        ggplot2::theme_minimal()

    result <- list(
        plots = plots,
        r_squared = r_squared,
        fitted = fitted_vals,
        residuals = residuals
    )

    return(result)
}

# =============================================================================
# Distribution-Specific Functions
# =============================================================================

#' PPC for Bernoulli imputation model
#'
#' @param fit CmdStanMCMC fit object
#' @param data Data frame
#' @param imputation_covariates Character vector of imputation model covariates
#' @param n_draws Number of posterior draws
#' @return List with calibration plot, ROC metrics
#' @export
ppc_imputation_bernoulli <- function(fit, data, imputation_covariates, n_draws = 1000) {
  ppc_data <- get_ppc_data(data, imputation_covariates)
  x_obs <- ppc_data$x
  Z_obs <- ppc_data$Z
  n_obs <- length(x_obs)

  draws <- posterior::as_draws_df(fit$draws())
  draw_idx <- sample(nrow(draws), min(n_draws, nrow(draws)))

  alpha_imp <- draws$alpha_imputation[draw_idx]
  gamma_mat <- extract_gamma_matrix(draws, draw_idx)

  prob_mat <- matrix(NA, nrow = length(draw_idx), ncol = n_obs)
  for (m in seq_along(draw_idx)) {
    eta <- alpha_imp[m]
    if (ncol(Z_obs) > 0) {
      eta <- eta + as.vector(Z_obs %*% gamma_mat[m, ])
    }
    prob_mat[m, ] <- stats::plogis(eta)
  }
  pred_prob <- colMeans(prob_mat)

  n_bins <- 10
  breaks <- seq(0, 1, length.out = n_bins + 1)

  calib_df <- data.frame(x = x_obs, prob = pred_prob) |>
    dplyr::mutate(
      bin = cut(prob, breaks = breaks, include.lowest = TRUE)
    ) |>
    dplyr::summarize(
      n = dplyr::n(),
      observed = mean(x),
      predicted = mean(prob),
      se = sqrt(observed * (1 - observed) / n),
      .by = bin
    ) |>
    dplyr::filter(!is.na(bin))

  p_calib <- ggplot2::ggplot(
    calib_df,
    ggplot2::aes(x = predicted, y = observed)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = pmax(0, observed - 1.96 * se), ymax = pmin(1, observed + 1.96 * se)),
      width = 0.02, color = "steelblue"
    ) +
    ggplot2::geom_point(ggplot2::aes(size = n), color = "steelblue") +
    ggplot2::coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = "Calibration: Bernoulli Imputation Model",
      subtitle = "Points on diagonal = well calibrated",
      x = "Mean Predicted P(x=1)",
      y = "Observed Proportion",
      size = "N"
    ) +
    ggplot2::theme_minimal()

  brier <- mean((pred_prob - x_obs)^2)

  thresholds <- seq(0, 1, by = 0.01)
  roc_data <- lapply(thresholds, \(thresh) {
    pred <- as.integer(pred_prob >= thresh)
    tp <- sum(pred == 1L & x_obs == 1L)
    fp <- sum(pred == 1L & x_obs == 0L)
    tn <- sum(pred == 0L & x_obs == 0L)
    fn <- sum(pred == 0L & x_obs == 1L)

    tpr <- if ((tp + fn) > 0L) tp / (tp + fn) else 0
    fpr <- if ((fp + tn) > 0L) fp / (fp + tn) else 0

    data.frame(threshold = thresh, tpr = tpr, fpr = fpr)
  }) |>
    dplyr::bind_rows()

  roc_data <- roc_data |> dplyr::arrange(fpr, tpr)
  auc <- sum(diff(roc_data$fpr) * (roc_data$tpr[-1] + roc_data$tpr[-nrow(roc_data)]) / 2)
  auc <- abs(auc)

  p_roc <- ggplot2::ggplot(roc_data, ggplot2::aes(x = fpr, y = tpr)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = "ROC Curve: Bernoulli Imputation Model",
      subtitle = sprintf("AUC = %.3f", auc),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    ggplot2::theme_minimal()

  return(list(
    calibration_plot = p_calib,
    roc_plot = p_roc,
    auc = auc,
    brier_score = brier,
    calibration_data = calib_df
  ))
}

#' PPC for count imputation model (negative binomial)
#'
#' @param fit CmdStanMCMC fit object
#' @param data Data frame
#' @param imputation_covariates Character vector of imputation model covariates
#' @param n_draws Number of posterior draws
#' @return List with rootogram and dispersion check
#' @export
ppc_imputation_count <- function(fit, data, imputation_covariates, n_draws = 500) {
  ppc_data <- get_ppc_data(data, imputation_covariates)
  x_obs <- ppc_data$x
  Z_obs <- ppc_data$Z
  n_obs <- length(x_obs)

  draws <- posterior::as_draws_df(fit$draws())
  draw_idx <- sample(nrow(draws), min(n_draws, nrow(draws)))

  alpha_imp <- draws$alpha_imputation[draw_idx]
  phi_col <- intersect(c("phi", "phi_imputation"), names(draws))
  if (length(phi_col) == 0L) {
    cli::cli_abort("Dispersion parameter {.var phi} not found in posterior draws.")
  }
  phi <- draws[[phi_col[1L]]][draw_idx]
  gamma_mat <- extract_gamma_matrix(draws, draw_idx)

  x_rep_list <- vector("list", length(draw_idx))
  for (m in seq_along(draw_idx)) {
    eta <- alpha_imp[m]
    if (ncol(Z_obs) > 0) {
      eta <- eta + as.vector(Z_obs %*% gamma_mat[m, ])
    }
    mu <- exp(eta)
    x_rep_list[[m]] <- stats::rnbinom(n_obs, mu = mu, size = phi[m])
  }
  x_rep_all <- unlist(x_rep_list)

  max_count <- min(max(c(x_obs, x_rep_all), na.rm = TRUE), 50)
  count_range <- 0:max_count

  obs_freq <- as.numeric(table(factor(x_obs, levels = count_range))) / n_obs
  exp_freq <- as.numeric(table(factor(x_rep_all, levels = count_range))) / length(x_rep_all)

  rootogram_df <- data.frame(
    count = count_range,
    observed = obs_freq,
    expected = exp_freq
  )

  p_rootogram <- ggplot2::ggplot(rootogram_df, ggplot2::aes(x = count)) +
    ggplot2::geom_col(
      ggplot2::aes(y = sqrt(observed)),
      fill = "steelblue", width = 0.7, alpha = 0.7
    ) +
    ggplot2::geom_point(ggplot2::aes(y = sqrt(expected)), color = "red", size = 2) +
    ggplot2::geom_line(ggplot2::aes(y = sqrt(expected)), color = "red") +
    ggplot2::labs(
      title = "Rootogram: Count Imputation Model",
      subtitle = "Bars = observed; Red = expected (sqrt scale)",
      x = "Count",
      y = expression(sqrt("Frequency"))
    ) +
    ggplot2::theme_minimal()

  obs_mean <- mean(x_obs)
  obs_var <- stats::var(x_obs)
  obs_dispersion <- obs_var / obs_mean

  rep_dispersions <- sapply(x_rep_list, \(x) stats::var(x) / mean(x))

  p_dispersion <- ggplot2::ggplot(
    data.frame(d = rep_dispersions),
    ggplot2::aes(x = d)
  ) +
    ggplot2::geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = obs_dispersion, color = "red", linewidth = 1.5) +
    ggplot2::labs(
      title = "Dispersion Check",
      subtitle = sprintf(
        "Observed = %.2f, Predicted median = %.2f",
        obs_dispersion, stats::median(rep_dispersions)
      ),
      x = "Variance / Mean",
      y = "Count"
    ) +
    ggplot2::theme_minimal()

  disp_p <- 2 * min(
    mean(rep_dispersions <= obs_dispersion),
    mean(rep_dispersions >= obs_dispersion)
  )

  return(list(
    rootogram = p_rootogram,
    dispersion_plot = p_dispersion,
    observed_mean = obs_mean,
    observed_dispersion = obs_dispersion,
    predicted_dispersion = stats::median(rep_dispersions),
    dispersion_p_value = disp_p
  ))
}

#' PPC for beta-binomial imputation model
#'
#' @param fit CmdStanMCMC fit object
#' @param data Data frame
#' @param imputation_covariates Character vector of imputation model covariates
#' @param n_draws Number of posterior draws
#' @return List with diagnostics
#' @export
ppc_imputation_beta_binomial <- function(fit, data, imputation_covariates, n_draws = 500) {
  ppc_data <- get_ppc_data(data, imputation_covariates)
  x_obs <- ppc_data$x
  Z_obs <- ppc_data$Z
  n_obs <- length(x_obs)

  n_trials <- data |>
    dplyr::summarize(n_trials = dplyr::first(n_trials), .by = id) |>
    dplyr::filter(id %in% ppc_data$id) |>
    dplyr::arrange(match(id, ppc_data$id)) |>
    dplyr::pull(n_trials)

  draws <- posterior::as_draws_df(fit$draws())
  draw_idx <- sample(nrow(draws), min(n_draws, nrow(draws)))

  alpha_imp <- draws$alpha_imputation[draw_idx]
  phi_col <- intersect(c("phi", "phi_imputation"), names(draws))
  if (length(phi_col) == 0L) {
    cli::cli_abort("Dispersion parameter {.var phi} not found in posterior draws.")
  }
  phi <- draws[[phi_col[1L]]][draw_idx]
  gamma_mat <- extract_gamma_matrix(draws, draw_idx)

  x_rep_list <- vector("list", length(draw_idx))
  for (m in seq_along(draw_idx)) {
    eta <- alpha_imp[m]
    if (ncol(Z_obs) > 0) {
      eta <- eta + as.vector(Z_obs %*% gamma_mat[m, ])
    }
    mu <- stats::plogis(eta)

    alpha_bb <- mu * phi[m]
    beta_bb <- (1 - mu) * phi[m]

    p_draw <- stats::rbeta(n_obs, alpha_bb, beta_bb)
    x_rep_list[[m]] <- stats::rbinom(n_obs, size = n_trials, prob = p_draw)
  }
  x_rep_all <- unlist(x_rep_list)

  max_count <- max(c(x_obs, n_trials))
  count_range <- 0:max_count

  obs_freq <- as.numeric(table(factor(x_obs, levels = count_range))) / n_obs
  exp_freq <- as.numeric(table(factor(x_rep_all, levels = count_range))) / length(x_rep_all)

  rootogram_df <- data.frame(
    count = count_range,
    observed = obs_freq,
    expected = exp_freq
  )

  p_rootogram <- ggplot2::ggplot(rootogram_df, ggplot2::aes(x = count)) +
    ggplot2::geom_col(
      ggplot2::aes(y = sqrt(observed)),
      fill = "steelblue", width = 0.7, alpha = 0.7
    ) +
    ggplot2::geom_point(ggplot2::aes(y = sqrt(expected)), color = "red", size = 2) +
    ggplot2::geom_line(ggplot2::aes(y = sqrt(expected)), color = "red") +
    ggplot2::labs(
      title = "Rootogram: Beta-Binomial Imputation Model",
      subtitle = "Bars = observed; Red = expected (sqrt scale)",
      x = "Count",
      y = expression(sqrt("Frequency"))
    ) +
    ggplot2::theme_minimal()

  prop_obs <- x_obs / n_trials
  prop_rep <- lapply(x_rep_list, \(x) x / n_trials)

  obs_prop_var <- stats::var(prop_obs)
  rep_prop_vars <- sapply(prop_rep, stats::var)

  var_p <- 2 * min(
    mean(rep_prop_vars <= obs_prop_var),
    mean(rep_prop_vars >= obs_prop_var)
  )

  return(list(
    rootogram = p_rootogram,
    observed_proportion_mean = mean(prop_obs),
    observed_proportion_var = obs_prop_var,
    predicted_proportion_var = stats::median(rep_prop_vars),
    variance_p_value = var_p
  ))
}

# =============================================================================
# Main Wrapper Function
# =============================================================================

#' Comprehensive posterior predictive checks for imputation model
#'
#' @param fit CmdStanMCMC fit object from fit_stan_model()
#' @param data Data frame used for fitting
#' @param imputation_covariates Character vector of imputation model covariate names
#' @param distribution Imputation distribution used
#' @param n_draws Number of posterior draws for simulations
#' @return List with all diagnostic results and plots
#' @export
ppc_imputation <- function(
  fit,
  data,
  imputation_covariates,
  distribution = c("normal", "bernoulli", "negative_binomial", "beta_binomial"),
  n_draws = 1000
) {
  distribution <- match.arg(distribution)

  cli::cli_h1("Imputation Model Diagnostics")
  cli::cli_alert_info("Distribution: {distribution}")

  ppc_data <- get_ppc_data(data, imputation_covariates)
  observed <- data.frame(id = ppc_data$id, x = ppc_data$x) # reconstruct for compatibility if needed or just use ppc_data
  cli::cli_alert_info("Observed x: n = {nrow(observed)}")

  results <- list()

  if (distribution == "normal") {
    cli::cli_h2("Density Comparison")
    results$density <- ppc_imputation_density(fit, data, imputation_covariates, n_draws = 100)

    cli::cli_h2("Summary Statistics")
    results$mean <- ppc_imputation_stat(fit, data, imputation_covariates,
      stat = mean, stat_name = "Mean", n_draws = n_draws
    )
    results$sd <- ppc_imputation_stat(fit, data, imputation_covariates,
      stat = stats::sd, stat_name = "SD", n_draws = n_draws
    )
    results$median <- ppc_imputation_stat(fit, data, imputation_covariates,
      stat = stats::median, stat_name = "Median", n_draws = n_draws
    )

    cli::cli_h2("Prediction Intervals")
    results$intervals <- ppc_imputation_intervals(fit, data, imputation_covariates,
      prob = 0.95, n_draws = n_draws
    )

    cli::cli_h2("Residual Diagnostics")
    results$residuals <- ppc_imputation_residuals(fit, data, imputation_covariates)

    cli::cli_h2("Mean Model Diagnostics")
    results$mean_diagnostics <- ppc_imputation_mean_diagnostics(
      fit, data, imputation_covariates, n_bins = 10
    )
  } else if (distribution == "bernoulli") {
    cli::cli_h2("Calibration and Discrimination")
    results$bernoulli <- ppc_imputation_bernoulli(fit, data, imputation_covariates, n_draws = n_draws)

    results$mean <- list(
      observed = mean(observed$x),
      p_value = NA
    )
  } else if (distribution == "negative_binomial") {
    cli::cli_h2("Count Distribution Check")
    results$count <- ppc_imputation_count(fit, data, imputation_covariates, n_draws = n_draws)

    results$mean <- ppc_imputation_stat(fit, data, imputation_covariates,
      stat = mean, stat_name = "Mean", n_draws = n_draws
    )
  } else if (distribution == "beta_binomial") {
    cli::cli_h2("Beta-Binomial Distribution Check")
    results$beta_binomial <- ppc_imputation_beta_binomial(fit, data, imputation_covariates, n_draws = n_draws)
  }

  cli::cli_h1("Summary")

  issues <- character(0)

  if (distribution == "normal") {
    cli::cli_alert_info("Mean check: observed = {round(results$mean$observed, 3)}, p = {round(results$mean$p_value, 3)}")
    cli::cli_alert_info("SD check: observed = {round(results$sd$observed, 3)}, p = {round(results$sd$p_value, 3)}")
    cli::cli_alert_info("95% interval coverage: {round(100 * results$intervals$coverage, 1)}%")
    cli::cli_alert_info("Shapiro-Wilk p-value: {round(results$residuals$diagnostics$shapiro_p, 3)}")
    cli::cli_alert_info("Mean model R-squared: {round(results$mean_diagnostics$r_squared, 3)}")

    if (results$mean$p_value < 0.05) {
      issues <- c(issues, "Mean of observed x differs from predicted (p < 0.05)")
    }
    if (results$sd$p_value < 0.05) {
      issues <- c(issues, "SD of observed x differs from predicted (p < 0.05)")
    }
    if (results$intervals$coverage < 0.85) {
      issues <- c(issues, sprintf(
        "Low interval coverage: %.0f%% (expected ~95%%)",
        100 * results$intervals$coverage
      ))
    }
    if (results$residuals$diagnostics$shapiro_p < 0.01) {
      issues <- c(issues, "Residuals significantly non-normal (Shapiro p < 0.01)")
    }
  } else if (distribution == "bernoulli") {
    cli::cli_alert_info("Brier score: {round(results$bernoulli$brier_score, 3)}")
    cli::cli_alert_info("AUC: {round(results$bernoulli$auc, 3)}")

    if (results$bernoulli$auc < 0.6) {
      issues <- c(issues, "Low AUC suggests poor predictive ability of imputation model")
    }
  } else if (distribution == "negative_binomial") {
    cli::cli_alert_info("Mean check: p = {round(results$mean$p_value, 3)}")
    cli::cli_alert_info("Dispersion check: observed = {round(results$count$observed_dispersion, 2)}, predicted = {round(results$count$predicted_dispersion, 2)}, p = {round(results$count$dispersion_p_value, 3)}")

    if (results$mean$p_value < 0.05) {
      issues <- c(issues, "Mean of observed x differs from predicted (p < 0.05)")
    }
    if (results$count$dispersion_p_value < 0.05) {
      issues <- c(issues, "Dispersion of observed x differs from predicted (p < 0.05)")
    }
  } else if (distribution == "beta_binomial") {
    cli::cli_alert_info("Proportion Mean: observed = {round(results$beta_binomial$observed_proportion_mean, 3)}")
    cli::cli_alert_info("Variance check: observed = {round(results$beta_binomial$observed_proportion_var, 4)}, predicted = {round(results$beta_binomial$predicted_proportion_var, 4)}, p = {round(results$beta_binomial$variance_p_value, 3)}")

    if (results$beta_binomial$variance_p_value < 0.05) {
      issues <- c(issues, "Variance of observed proportions differs from predicted (p < 0.05)")
    }
  }

  if (length(issues) > 0) {
    cli::cli_h2("Potential Issues")
    for (issue in issues) {
      cli::cli_alert_warning(issue)
    }
    cli::cli_alert_info("Consider: different distribution, non-linear terms, or additional covariates")
  } else {
    cli::cli_alert_success("No obvious misspecification detected")
  }

  results$issues <- issues
  results$distribution <- distribution

  class(results) <- c("ppc_imputation", class(results))

  return(results)
}

#' Print method for ppc_imputation results
#' @param x ppc_imputation results
#' @param ... Additional arguments (unused)
#' @export
print.ppc_imputation <- function(x, ...) {
  cat("Imputation Model PPC Results\n")
  cat("Distribution:", x$distribution, "\n")
  cat("Issues found:", length(x$issues), "\n")
  if (length(x$issues) > 0) {
    for (issue in x$issues) {
      cat("  -", issue, "\n")
    }
  }
  invisible(x)
}
