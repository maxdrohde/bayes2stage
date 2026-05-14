################################################################################
# Relative Efficiency Heatmap
################################################################################

library(ggplot2)
source("simulation_config.R")

SELECTED_PARAMETER <- "beta_x"

# ------------------------------------------------------------------------------
# Load and compute relative efficiency
# ------------------------------------------------------------------------------

sim_results <- arrow::open_dataset("processed_data/combined_data") |>
    dplyr::filter(parameter == SELECTED_PARAMETER) |>
    dplyr::collect()

# Compute variance for each (setting, type), then RE relative to SRS (no imp)
variance_by_type <- sim_results |>
    dplyr::summarise(variance = var(estimate, na.rm = TRUE), .by = c(sim_setting, type))

baseline_var <- variance_by_type |>
    dplyr::filter(type == "SRS\n(no imp)") |>
    dplyr::select(sim_setting, baseline_var = variance)

plot_data <- get_simulation_grid() |>
    dplyr::mutate(sim_setting = dplyr::row_number()) |>
    dplyr::select(sim_setting, beta_z, rand_intercept_sd, gamma1) |>
    dplyr::left_join(variance_by_type, by = "sim_setting") |>
    dplyr::filter(type != "SRS\n(no imp)") |>
    dplyr::left_join(baseline_var, by = "sim_setting") |>
    dplyr::mutate(
        rel_eff = baseline_var / variance,
        log2_re = log2(rel_eff),
        type = factor(type, levels = c("SRS", "ODS", "BDS", "ACML ODS")),
        beta_z_label = factor(
            paste0("\u03b2z = ", beta_z),
            levels = paste0("\u03b2z = ", c(0, -2, -4))
        )
    )

# ------------------------------------------------------------------------------
# Plot
# ------------------------------------------------------------------------------

p <- ggplot(plot_data, aes(x = factor(gamma1), y = factor(rand_intercept_sd))) +
    geom_tile(aes(fill = log2_re), color = "white", linewidth = 0.5) +
    geom_text(
        aes(label = sprintf("%.2f", rel_eff), color = abs(log2_re) > 0.7),
        size = 3.5, fontface = "bold"
    ) +
    scale_fill_gradientn(
        colours = rev(scico::scico(256, palette = "vik")),
        limits = c(-2, 2),
        oob = scales::squish,
        name = "Relative\nEfficiency",
        breaks = c(-2, -1, 0, 1, 2),
        labels = c("0.25", "0.5", "1", "2", "4")
    ) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
    facet_grid(rows = vars(beta_z_label), cols = vars(type)) +
    coord_equal() +
    labs(
        title = glue::glue("Relative Efficiency: {SELECTED_PARAMETER}"),
        subtitle = "Relative efficiency compared to SRS without imputation (values > 1 indicate improved efficiency)",
        x = expression(gamma[1] ~ "(z-x association)"),
        y = expression(sigma[b0] ~ "(random intercept SD)"),
        caption = paste0(
            "**Varied:** beta_z = {0, -2, -4} (rows) | sigma_b0 = {1, 3, 6, 9} (y-axis) | gamma1 = {0, 0.5, 1, 2} (x-axis) | 48 settings total<br>",
            "**Fixed:** N = 300 | 25% stage-2 sampling | 4-6 visits/subject | ODS/BDS: 10th/90th percentile cutoffs, 40%/20%/40% allocation<br>",
            "**Outcome:** y = 75 - 0.5x + beta_z(z) - t - 0.5(xt) + b0 + b1(t) + e | sigma_b1 = 1.25 | rho = 0 | sigma_eps = 1<br>",
            "**DGP:** x ~ Bernoulli(0.25) | z | x ~ N(0.25 + gamma1(x), 1)"
        )
    ) +
    cowplot::theme_cowplot(font_size = 11, font_family = "Source Sans Pro") +
    theme(
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(face = "bold", size = 10),
        legend.position = "right",
        legend.key.height = unit(2, "cm"),
        panel.grid = element_blank(),
        panel.spacing = unit(0.3, "lines"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        plot.caption = ggtext::element_markdown(size = 10, hjust = 0, lineheight = 1.4),
        plot.caption.position = "plot",
        plot.margin = margin(10, 10, 15, 10)
    )

# ------------------------------------------------------------------------------
# Save
# ------------------------------------------------------------------------------

fs::dir_create("plots")
ggsave(
    glue::glue("plots/relative_efficiency_heatmap_{SELECTED_PARAMETER}.png"),
    plot = p, width = 14, height = 11, dpi = 300, bg = "white", device = ragg::agg_png
)
cli::cli_alert_success("Saved heatmap to plots/")
