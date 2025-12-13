#' @title MCMC Diagnostic Plots
#' @name mcmc_diagnostics
#' @description Functions for diagnosing MCMC sampling quality
NULL

# =============================================================================
# Default Parameters
# =============================================================================

#' Default key parameters for diagnostics
#' @keywords internal
default_key_params <- function() {
    c("alpha_main", "beta_t", "beta_x", "beta_x_t_interaction",
      "beta[1]", "alpha_imputation", "gamma[1]",
      "sigma_main", "sigma_imputation", "sigma_re[1]", "sigma_re[2]")
}

# =============================================================================
# Individual Plot Functions
# =============================================================================

#' Trace plots for MCMC diagnostics
#'
#' @param fit CmdStanMCMC fit object from fit_stan_model()
#' @param parameters Character vector of parameter names to plot.
#'   If NULL, uses default key parameters.
#' @return A ggplot object
#' @export
plot_trace <- function(fit, parameters = NULL) {
    if (is.null(parameters)) {
        parameters <- default_key_params()
    }

    draws <- posterior::as_draws_array(fit$draws())

    available <- posterior::variables(draws)
    parameters <- intersect(parameters, available)

    if (length(parameters) == 0) {
        cli::cli_abort("No matching parameters found in draws")
    }

    bayesplot::mcmc_trace(draws, pars = parameters) +
        ggplot2::labs(title = "Trace Plots")
}

#' Pairs plots for diagnosing posterior geometry
#'
#' @param fit CmdStanMCMC fit object
#' @param parameters Character vector of parameter names to plot.
#'   If NULL, uses a subset of key parameters.
#' @return A ggplot object
#' @export
plot_pairs <- function(fit, parameters = NULL) {
    if (is.null(parameters)) {
        parameters <- c("beta_x", "sigma_re[1]", "sigma_re[2]", "sigma_main")
    }

    draws <- posterior::as_draws_array(fit$draws())

    available <- posterior::variables(draws)
    parameters <- intersect(parameters, available)

    if (length(parameters) < 2) {
        cli::cli_abort("Need at least 2 matching parameters for pairs plot")
    }

    np <- bayesplot::nuts_params(fit)

    bayesplot::mcmc_pairs(
        draws,
        pars = parameters,
        off_diag_args = list(size = 0.5, alpha = 0.3),
        np = np
    ) +
        ggplot2::labs(title = "Pairs Plot (divergences highlighted)")
}

#' Rank histogram plots for convergence diagnostics
#'
#' @param fit CmdStanMCMC fit object
#' @param parameters Character vector of parameter names to plot.
#'   If NULL, uses default key parameters.
#' @return A ggplot object
#' @export
plot_rank <- function(fit, parameters = NULL) {
    if (is.null(parameters)) {
        parameters <- default_key_params()
    }

    draws <- posterior::as_draws_array(fit$draws())

    available <- posterior::variables(draws)
    parameters <- intersect(parameters, available)

    if (length(parameters) == 0) {
        cli::cli_abort("No matching parameters found in draws")
    }

    bayesplot::mcmc_rank_histogram(draws, pars = parameters) +
        ggplot2::labs(title = "Rank Histograms (uniform = good mixing)")
}

#' Density plots with chain overlay
#'
#' @param fit CmdStanMCMC fit object
#' @param parameters Character vector of parameter names to plot.
#'   If NULL, uses default key parameters.
#' @return A ggplot object
#' @export
plot_density <- function(fit, parameters = NULL) {
    if (is.null(parameters)) {
        parameters <- default_key_params()
    }

    draws <- posterior::as_draws_array(fit$draws())

    available <- posterior::variables(draws)
    parameters <- intersect(parameters, available)

    if (length(parameters) == 0) {
        cli::cli_abort("No matching parameters found in draws")
    }

    bayesplot::mcmc_dens_overlay(draws, pars = parameters) +
        ggplot2::labs(title = "Posterior Densities by Chain")
}

#' Posterior interval plots
#'
#' @param fit CmdStanMCMC fit object
#' @param parameters Character vector of parameter names to plot.
#'   If NULL, uses default key parameters.
#' @param prob Inner probability interval (default: 0.5)
#' @param prob_outer Outer probability interval (default: 0.95)
#' @return A ggplot object
#' @export
plot_intervals <- function(fit, parameters = NULL, prob = 0.5, prob_outer = 0.95) {
    if (is.null(parameters)) {
        parameters <- default_key_params()
    }

    draws <- posterior::as_draws_array(fit$draws())

    available <- posterior::variables(draws)
    parameters <- intersect(parameters, available)

    if (length(parameters) == 0) {
        cli::cli_abort("No matching parameters found in draws")
    }

    bayesplot::mcmc_intervals(draws, pars = parameters,
                              prob = prob, prob_outer = prob_outer) +
        ggplot2::labs(title = "Posterior Intervals")
}

#' Energy diagnostic plot (NUTS-specific)
#'
#' @param fit CmdStanMCMC fit object
#' @return A ggplot object
#' @export
plot_energy <- function(fit) {
    np <- bayesplot::nuts_params(fit)

    bayesplot::mcmc_nuts_energy(np) +
        ggplot2::labs(title = "NUTS Energy Diagnostic",
                      subtitle = "Overlapping distributions = good")
}

#' Pareto k diagnostic plot for divergences
#'
#' @param fit CmdStanMCMC fit object
#' @param parameters Character vector of parameter names.
#'   If NULL, uses parameters related to random effects.
#' @return A ggplot object
#' @export
plot_pareto_k <- function(fit, parameters = NULL) {
    if (is.null(parameters)) {
        parameters <- c("sigma_re[1]", "sigma_re[2]")
    }

    draws <- posterior::as_draws_array(fit$draws())
    np <- bayesplot::nuts_params(fit)

    available <- posterior::variables(draws)
    parameters <- intersect(parameters, available)

    if (length(parameters) == 0) {
        cli::cli_abort("No matching parameters found in draws")
    }

    bayesplot::mcmc_parcoord(draws, pars = parameters, np = np) +
        ggplot2::labs(title = "Parallel Coordinates (divergences in red)")
}

# =============================================================================
# Comprehensive Diagnostics
# =============================================================================

#' Comprehensive MCMC diagnostics
#'
#' Runs a full suite of MCMC diagnostics and returns plots and numeric summaries.
#'
#' @param fit CmdStanMCMC fit object from fit_stan_model()
#' @param parameters Character vector of parameter names to diagnose.
#'   If NULL, uses default key parameters.
#' @param print_summary Logical; print diagnostic summary to console (default: TRUE)
#' @return A list with plots and diagnostic summaries
#' @export
diagnose_fit <- function(fit, parameters = NULL, print_summary = TRUE) {
    if (is.null(parameters)) {
        parameters <- default_key_params()
    }

    draws <- posterior::as_draws_array(fit$draws())
    available <- posterior::variables(draws)
    parameters <- intersect(parameters, available)

    if (length(parameters) == 0) {
        cli::cli_abort("No matching parameters found in draws")
    }

    diag <- fit$diagnostic_summary()
    summ <- fit$summary(variables = parameters)

    if (print_summary) {
        cli::cli_h1("MCMC Diagnostics")

        cli::cli_h2("Sampling Diagnostics")
        cli::cli_alert_info("Divergences: {sum(diag$num_divergent)}")
        cli::cli_alert_info("Max treedepth exceeded: {sum(diag$num_max_treedepth)}")
        cli::cli_alert_info("E-BFMI: {paste(round(diag$ebfmi, 3), collapse = ', ')}")

        cli::cli_h2("Parameter Diagnostics")
        cli::cli_alert_info("Min ESS (bulk): {round(min(summ$ess_bulk, na.rm = TRUE), 1)}")
        cli::cli_alert_info("Min ESS (tail): {round(min(summ$ess_tail, na.rm = TRUE), 1)}")
        cli::cli_alert_info("Max Rhat: {round(max(summ$rhat, na.rm = TRUE), 4)}")

        if (sum(diag$num_divergent) > 0) {
            cli::cli_alert_warning("Divergences detected! Check pairs plot for problematic parameters.")
        }
        if (max(summ$rhat, na.rm = TRUE) > 1.01) {
            cli::cli_alert_warning("High Rhat detected! Chains may not have converged.")
        }
        if (min(summ$ess_bulk, na.rm = TRUE) < 400) {
            cli::cli_alert_warning("Low ESS detected! Consider more iterations.")
        }
        if (any(diag$ebfmi < 0.3)) {
            cli::cli_alert_warning("Low E-BFMI detected! May indicate poor posterior geometry.")
        }
    }

    plots <- list()

    cli::cli_progress_step("Generating trace plots...")
    plots$trace <- plot_trace(fit, parameters)

    cli::cli_progress_step("Generating rank histograms...")
    plots$rank <- plot_rank(fit, parameters)

    cli::cli_progress_step("Generating density plots...")
    plots$density <- plot_density(fit, parameters)

    cli::cli_progress_step("Generating interval plot...")
    plots$intervals <- plot_intervals(fit, parameters)

    cli::cli_progress_step("Generating energy plot...")
    plots$energy <- plot_energy(fit)

    pairs_params <- intersect(c("beta_x", "sigma_re[1]", "sigma_re[2]", "sigma_main"), available)
    if (length(pairs_params) >= 2) {
        cli::cli_progress_step("Generating pairs plot...")
        plots$pairs <- plot_pairs(fit, pairs_params)
    }

    cli::cli_progress_done()

    results <- list(
        plots = plots,
        summary = summ,
        diagnostics = list(
            divergences = sum(diag$num_divergent),
            max_treedepth = sum(diag$num_max_treedepth),
            ebfmi = diag$ebfmi,
            min_ess_bulk = min(summ$ess_bulk, na.rm = TRUE),
            min_ess_tail = min(summ$ess_tail, na.rm = TRUE),
            max_rhat = max(summ$rhat, na.rm = TRUE)
        ),
        warnings = character(0)
    )

    if (sum(diag$num_divergent) > 0) {
        results$warnings <- c(results$warnings, "Divergences detected")
    }
    if (max(summ$rhat, na.rm = TRUE) > 1.01) {
        results$warnings <- c(results$warnings, "High Rhat (> 1.01)")
    }
    if (min(summ$ess_bulk, na.rm = TRUE) < 400) {
        results$warnings <- c(results$warnings, "Low ESS (< 400)")
    }
    if (any(diag$ebfmi < 0.3)) {
        results$warnings <- c(results$warnings, "Low E-BFMI (< 0.3)")
    }

    class(results) <- c("bayes2stage_diagnostics", class(results))

    return(results)
}

#' Print method for bayes2stage diagnostics
#' @param x bayes2stage_diagnostics object
#' @param ... Additional arguments (unused)
#' @export
print.bayes2stage_diagnostics <- function(x, ...) {
    cat("MCMC Diagnostics Summary\n")
    cat(strrep("-", 40), "\n")
    cat("Divergences:", x$diagnostics$divergences, "\n")
    cat("Max treedepth:", x$diagnostics$max_treedepth, "\n")
    cat("Min ESS (bulk):", round(x$diagnostics$min_ess_bulk, 1), "\n")
    cat("Max Rhat:", round(x$diagnostics$max_rhat, 4), "\n")

    if (length(x$warnings) > 0) {
        cat("\nWarnings:\n")
        for (w in x$warnings) {
            cat("  -", w, "\n")
        }
    } else {
        cat("\nNo warnings.\n")
    }

    cat("\nAvailable plots: ", paste(names(x$plots), collapse = ", "), "\n")
    invisible(x)
}

#' Quick funnel diagnostic for hierarchical parameters
#'
#' Creates a scatter plot of log(sigma) vs random effects to diagnose funnel geometry.
#'
#' @param fit CmdStanMCMC fit object
#' @param sigma_param Name of the sigma parameter (default: `"sigma_re[1]"`)
#' @param effect_param Pattern for random effect parameters (default: `"re\\[1,"`)
#' @param n_effects Number of random effects to plot (default: 5)
#' @return A ggplot object
#' @export
plot_funnel <- function(fit, sigma_param = "sigma_re[1]",
                        effect_param = "re\\[1,", n_effects = 5) {
    draws_df <- posterior::as_draws_df(fit$draws())

    if (!sigma_param %in% names(draws_df)) {
        cli::cli_abort("Parameter {sigma_param} not found in draws")
    }

    effect_cols <- grep(effect_param, names(draws_df), value = TRUE)
    if (length(effect_cols) == 0) {
        cli::cli_abort("No parameters matching pattern {effect_param} found")
    }

    effect_cols <- effect_cols[seq_len(min(n_effects, length(effect_cols)))]

    plot_data <- lapply(effect_cols, function(col) {
        data.frame(
            log_sigma = log(draws_df[[sigma_param]]),
            effect = draws_df[[col]],
            param = col
        )
    }) |>
        dplyr::bind_rows()

    ggplot2::ggplot(plot_data, ggplot2::aes(x = log_sigma, y = effect)) +
        ggplot2::geom_point(alpha = 0.1, size = 0.5) +
        ggplot2::facet_wrap(~param, scales = "free_y") +
        ggplot2::labs(
            title = "Funnel Diagnostic",
            subtitle = "Narrow funnel at low sigma = potential sampling issues",
            x = paste0("log(", sigma_param, ")"),
            y = "Random Effect"
        ) +
        ggplot2::theme_minimal()
}
