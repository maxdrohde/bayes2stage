# Analysis Code for Testing Simulation Hypotheses
# bayes2stage package
# Generated: 2024-12-17

# ==============================================================================
# Setup
# ==============================================================================

library(ggplot2)
library(patchwork)

theme_set(
    cowplot::theme_cowplot(font_size = 12, font_family = "Source Sans Pro")
)

# ==============================================================================
# H1: Efficiency Gains Scale with Effect Size
# ==============================================================================

analyze_h1 <- function(results, true_values) {
    # Compute summary statistics by method and beta_z
    h1_summary <- results |>
        dplyr::filter(parameter == "beta_z") |>
        dplyr::left_join(true_values, by = c("sim_setting", "parameter")) |>
        dplyr::summarise(
            bias = mean(estimate) - unique(true_value),
            empirical_se = sd(estimate),
            n_sims = dplyr::n(),
            .by = c(type, true_value)
        )

    # Compute efficiency ratios
    h1_efficiency <- h1_summary |>
        dplyr::select(type, true_value, empirical_se) |>
        tidyr::pivot_wider(names_from = type, values_from = empirical_se) |>
        dplyr::mutate(
            abs_beta_z = abs(true_value)
        )

    # Plot: SE by method across beta_z
    p1 <- ggplot(h1_summary, aes(x = true_value, y = empirical_se, color = type)) +
        geom_line(linewidth = 1) +
        geom_point(size = 3) +
        scale_color_brewer(palette = "Set1") +
        labs(
            title = "Standard Error vs Effect Size",
            x = expression(beta[z]~"(true value)"),
            y = "Empirical SE",
            color = "Method"
        )

    result <- list(
        summary = h1_summary,
        efficiency = h1_efficiency,
        plot = p1
    )
    return(result)
}

# ==============================================================================
# H2: Asymmetry Around Zero
# ==============================================================================

analyze_h2 <- function(results, true_values) {
    h2_summary <- results |>
        dplyr::filter(parameter == "beta_z") |>
        dplyr::left_join(true_values, by = c("sim_setting", "parameter")) |>
        dplyr::filter(true_value != 0) |>
        dplyr::summarise(
            bias = mean(estimate) - unique(true_value),
            empirical_se = sd(estimate),
            .by = c(type, true_value)
        ) |>
        dplyr::mutate(
            abs_beta_z = abs(true_value),
            sign_beta_z = sign(true_value)
        )

    # Test for asymmetry
    h2_asymmetry <- h2_summary |>
        tidyr::pivot_wider(
            names_from = sign_beta_z,
            values_from = c(bias, empirical_se),
            names_sep = "_sign_"
        ) |>
        dplyr::mutate(
            bias_diff = `bias_sign_1` - `bias_sign_-1`,
            se_diff = `empirical_se_sign_1` - `empirical_se_sign_-1`
        )

    # Plot: Bias symmetry
    p1 <- ggplot(h2_summary, aes(x = true_value, y = bias, color = type)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        geom_line(linewidth = 1) +
        geom_point(size = 3) +
        scale_color_brewer(palette = "Set1") +
        labs(
            title = "Bias Symmetry Around Zero",
            x = expression(beta[z]~"(true value)"),
            y = "Bias",
            color = "Method"
        )

    result <- list(
        summary = h2_summary,
        asymmetry = h2_asymmetry,
        plot = p1
    )
    return(result)
}

# ==============================================================================
# H3: Null Effect Behavior
# ==============================================================================

analyze_h3 <- function(results, true_values) {
    h3_summary <- results |>
        dplyr::filter(parameter == "beta_z") |>
        dplyr::left_join(true_values, by = c("sim_setting", "parameter")) |>
        dplyr::filter(true_value == 0) |>
        dplyr::summarise(
            mean_estimate = mean(estimate),
            bias = mean(estimate),
            empirical_se = sd(estimate),
            model_se = mean(se),
            n_sims = dplyr::n(),
            .by = type
        )

    # Plot: Distribution of estimates at null
    p1 <- results |>
        dplyr::filter(parameter == "beta_z") |>
        dplyr::left_join(true_values, by = c("sim_setting", "parameter")) |>
        dplyr::filter(true_value == 0) |>
        ggplot(aes(x = estimate, fill = type)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
        geom_density(alpha = 0.5) +
        scale_fill_brewer(palette = "Set1") +
        labs(
            title = "Distribution of Estimates at Null Effect",
            x = expression(hat(beta)[z]),
            y = "Density",
            fill = "Method"
        )

    result <- list(
        summary = h3_summary,
        plot = p1
    )
    return(result)
}

# ==============================================================================
# H4: Bias-Variance Tradeoff
# ==============================================================================

analyze_h4 <- function(results, true_values) {
    h4_summary <- results |>
        dplyr::filter(parameter == "beta_z") |>
        dplyr::left_join(true_values, by = c("sim_setting", "parameter")) |>
        dplyr::summarise(
            bias = mean(estimate) - unique(true_value),
            variance = var(estimate),
            empirical_se = sd(estimate),
            mse = mean((estimate - unique(true_value))^2),
            .by = c(type, true_value)
        ) |>
        dplyr::mutate(
            abs_bias = abs(bias),
            mse_bias_component = bias^2,
            mse_variance_component = variance,
            pct_mse_from_bias = 100 * mse_bias_component / mse
        )

    # Plot: MSE decomposition
    p1 <- h4_summary |>
        dplyr::select(type, true_value, mse_bias_component, mse_variance_component) |>
        tidyr::pivot_longer(
            cols = dplyr::starts_with("mse_"),
            names_to = "component",
            values_to = "value",
            names_prefix = "mse_"
        ) |>
        ggplot(aes(x = true_value, y = value, fill = component)) +
        geom_col(position = "stack") +
        facet_wrap(~type) +
        scale_fill_manual(
            values = c("bias_component" = "#E41A1C", "variance_component" = "#377EB8"),
            labels = c("Bias^2", "Variance")
        ) +
        labs(
            title = "MSE Decomposition",
            x = expression(beta[z]),
            y = "MSE Component",
            fill = "Component"
        )

    result <- list(
        summary = h4_summary,
        plot = p1
    )
    return(result)
}

# ==============================================================================
# H5: Imputation Model Strength (gamma1 variation)
# ==============================================================================

analyze_h5 <- function(results, true_values, grid) {
    # Join with grid to get gamma1
    grid_gamma <- grid |>
        dplyr::mutate(sim_setting = dplyr::row_number()) |>
        dplyr::select(sim_setting, gamma1)

    h5_summary <- results |>
        dplyr::filter(parameter == "beta_z") |>
        dplyr::left_join(grid_gamma, by = "sim_setting") |>
        dplyr::summarise(
            empirical_se = sd(estimate),
            .by = c(type, gamma1)
        )

    # Plot: SE vs gamma1
    p1 <- ggplot(h5_summary, aes(x = gamma1, y = empirical_se, color = type)) +
        geom_line(linewidth = 1) +
        geom_point(size = 3) +
        scale_color_brewer(palette = "Set1") +
        labs(
            title = "Effect of Imputation Model Strength",
            x = expression(gamma[1]~"(z -> x association)"),
            y = "Empirical SE",
            color = "Method"
        )

    result <- list(
        summary = h5_summary,
        plot = p1
    )
    return(result)
}

# ==============================================================================
# Comprehensive Analysis Function
# ==============================================================================

analyze_all_hypotheses <- function(results, true_values, grid) {
    cli::cli_h1("Analyzing Simulation Hypotheses")

    output <- list()

    cli::cli_alert_info("H1: Efficiency scales with effect size")
    output$h1 <- analyze_h1(results, true_values)

    cli::cli_alert_info("H2: Asymmetry around zero")
    output$h2 <- analyze_h2(results, true_values)

    if (0 %in% true_values$true_value) {
        cli::cli_alert_info("H3: Null effect behavior")
        output$h3 <- analyze_h3(results, true_values)
    }

    cli::cli_alert_info("H4: Bias-variance tradeoff")
    output$h4 <- analyze_h4(results, true_values)

    if ("gamma1" %in% names(grid) && length(unique(grid$gamma1)) > 1) {
        cli::cli_alert_info("H5: Imputation model strength")
        output$h5 <- analyze_h5(results, true_values, grid)
    }

    cli::cli_alert_success("Analysis complete")
    return(output)
}

# ==============================================================================
# Summary Statistics Helper
# ==============================================================================

compute_simulation_metrics <- function(results, true_values, param_name) {
    results |>
        dplyr::filter(parameter == param_name) |>
        dplyr::left_join(true_values, by = c("sim_setting", "parameter")) |>
        dplyr::summarise(
            mean_estimate = mean(estimate),
            bias = mean(estimate) - unique(true_value),
            pct_bias = dplyr::if_else(
                unique(true_value) != 0,
                100 * (mean(estimate) - unique(true_value)) / abs(unique(true_value)),
                NA_real_
            ),
            empirical_se = sd(estimate),
            model_se = mean(se),
            se_bias_pct = 100 * (mean(se) - sd(estimate)) / sd(estimate),
            rmse = sqrt(mean((estimate - unique(true_value))^2)),
            mse = mean((estimate - unique(true_value))^2),
            n_sims = dplyr::n(),
            .by = c(type, true_value)
        )
}

# ==============================================================================
# Efficiency Table
# ==============================================================================

create_efficiency_table <- function(summary_data, reference_method = "ODS") {
    ref_se <- summary_data |>
        dplyr::filter(type == reference_method) |>
        dplyr::pull(empirical_se)

    summary_data |>
        dplyr::mutate(
            relative_efficiency = ref_se / empirical_se,
            pct_se_change = 100 * (empirical_se - ref_se) / empirical_se
        ) |>
        dplyr::arrange(empirical_se)
}

# ==============================================================================
# Publication-Ready Table
# ==============================================================================

format_results_table <- function(summary_data) {
    summary_data |>
        dplyr::transmute(
            Method = type,
            Bias = sprintf("%.3f", bias),
            `% Bias` = sprintf("%.1f%%", pct_bias),
            `Emp. SE` = sprintf("%.3f", empirical_se),
            `Model SE` = sprintf("%.3f", model_se),
            RMSE = sprintf("%.3f", rmse)
        )
}
