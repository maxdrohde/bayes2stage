################################################################################
# Relative Efficiency Table
#
# Creates a gt table showing relative efficiency vs SRS (no imputation)
# for all 48 simulation settings.
################################################################################

library(ggplot2)
library(gt)

source("simulation_config.R")

# ==============================================================================
# Configuration
# ==============================================================================

# Toggle: which parameter to analyze
SELECTED_PARAMETER <- "beta_x"

# ==============================================================================
# Load and Process Data
# ==============================================================================

cached_dataset <- "processed_data/combined_data"

if (!fs::dir_exists(cached_dataset)) {
    cli::cli_abort(
        "Data not found at {.path {cached_dataset}}. Run process_results.R first."
    )
}

sim_results <- arrow::open_dataset(cached_dataset) |>
    dplyr::filter(parameter == SELECTED_PARAMETER) |>
    dplyr::collect()

cli::cli_alert_success("Loaded {.val {nrow(sim_results)}} rows for parameter {.val {SELECTED_PARAMETER}}")

# Get simulation grid for setting labels
grid <- get_simulation_grid()

# ==============================================================================
# Compute Empirical Variance by Type and Setting
# ==============================================================================

# Compute empirical SE (SD of estimates) for each (setting, type) combination
empirical_se <- sim_results |>
    dplyr::summarise(
        se_emp = sd(estimate, na.rm = TRUE),
        n_sims = dplyr::n(),
        .by = c(sim_setting, type)
    ) |>
    dplyr::mutate(variance = se_emp^2)

# ==============================================================================
# Compute Relative Efficiency
# ==============================================================================

baseline_type <- "SRS\n(no imp)"

# Get baseline variance for each setting
baseline_var <- empirical_se |>
    dplyr::filter(type == baseline_type) |>
    dplyr::select(sim_setting, baseline_variance = variance)

# Compute relative efficiency for all other types
rel_eff <- empirical_se |>
    dplyr::filter(type != baseline_type) |>
    dplyr::left_join(baseline_var, by = "sim_setting") |>
    dplyr::mutate(
        rel_eff = baseline_variance / variance
    ) |>
    dplyr::select(sim_setting, type, rel_eff)

# Pivot to wide format: one column per type
rel_eff_wide <- rel_eff |>
    tidyr::pivot_wider(
        names_from = type,
        values_from = rel_eff
    )

# ==============================================================================
# Create Setting Labels
# ==============================================================================

setting_labels <- grid |>
    dplyr::mutate(
        sim_setting = dplyr::row_number()
    ) |>
    dplyr::select(sim_setting, beta_z, rand_intercept_sd, gamma1)

# Join setting labels with relative efficiency
table_data <- setting_labels |>
    dplyr::left_join(rel_eff_wide, by = "sim_setting") |>
    dplyr::mutate(
        # Create formatted row index
        `#` = sim_setting
    ) |>
    dplyr::select(
        `#`,
        `beta_z` = beta_z,
        `sigma_int` = rand_intercept_sd,
        `gamma1` = gamma1,
        SRS,
        ODS,
        BDS,
        `ACML ODS`
    )

# ==============================================================================
# Create gt Table
# ==============================================================================

# Define diverging color palette using RColorBrewer
# RdYlGn: red (inefficient) -> yellow (neutral) -> green (efficient)
rel_eff_cols <- c("SRS", "ODS", "BDS", "ACML ODS")

# Find range for color scaling
all_values <- unlist(table_data[, rel_eff_cols])
value_range <- range(all_values, na.rm = TRUE)

# Create symmetric range around 1 for balanced color scale
max_dev <- max(abs(value_range - 1))
color_min <- 1 - max_dev
color_max <- 1 + max_dev

setting_cols <- c("#", "beta_z", "sigma_int", "gamma1")

# Simulation settings description
settings_description <- paste0(
    "48 settings from full factorial: ",
    "&beta;<sub>z</sub> &isin; {-4, -2, 0}, ",
    "&sigma;<sub>int</sub> &isin; {1, 3, 6, 9}, ",
    "&gamma;<sub>1</sub> &isin; {0, 0.5, 1, 2}"
)

fixed_settings <- paste0(
    "Fixed: N=300, sampling fraction=0.25, ",
    "&beta;<sub>x</sub>=-0.5, &beta;<sub>t</sub>=-1, ",
    "&beta;<sub>x:t</sub>=-0.5, &sigma;<sub>slope</sub>=1.25"
)

gt_table <- table_data |>
    gt() |>
    tab_header(
        title = "Relative Efficiency vs SRS (no imputation)",
        subtitle = glue::glue("Parameter: {SELECTED_PARAMETER} | Values > 1 indicate higher efficiency")
    ) |>
    fmt_number(
        columns = dplyr::all_of(rel_eff_cols),
        decimals = 2
    ) |>
    data_color(
        columns = dplyr::all_of(rel_eff_cols),
        method = "numeric",
        palette = "RdYlGn",
        domain = c(color_min, color_max),
        na_color = "gray90"
    ) |>
    tab_spanner(
        label = "Setting",
        columns = dplyr::all_of(setting_cols)
    ) |>
    tab_spanner(
        label = "Relative Efficiency",
        columns = dplyr::all_of(rel_eff_cols)
    ) |>
    cols_align(
        align = "center",
        columns = dplyr::everything()
    ) |>
    cols_label(
        `#` = "#",
        beta_z = md("&beta;<sub>z</sub>"),
        sigma_int = md("&sigma;<sub>int</sub>"),
        gamma1 = md("&gamma;<sub>1</sub>")
    ) |>
    tab_options(
        table.font.size = px(12),
        heading.title.font.size = px(16),
        heading.subtitle.font.size = px(12),
        column_labels.font.weight = "bold",
        table_body.hlines.color = "gray90",
        table.border.top.color = "gray50",
        table.border.bottom.color = "gray50"
    ) |>
    tab_source_note(
        source_note = md(settings_description)
    ) |>
    tab_source_note(
        source_note = md(fixed_settings)
    ) |>
    tab_footnote(
        footnote = "Relative Efficiency = Var(SRS no imp) / Var(method)"
    )

# ==============================================================================
# Save Output
# ==============================================================================

output_dir <- "plots"
fs::dir_create(output_dir)

# Save as HTML
gtsave(
    gt_table,
    filename = glue::glue("{output_dir}/relative_efficiency_{SELECTED_PARAMETER}.html")
)
cli::cli_alert_success("Saved HTML table to {.path {output_dir}/relative_efficiency_{SELECTED_PARAMETER}.html}")

# Save as PNG
gtsave(
    gt_table,
    filename = glue::glue("{output_dir}/relative_efficiency_{SELECTED_PARAMETER}.png"),
    vwidth = 900,
    vheight = 1400
)
cli::cli_alert_success("Saved PNG table to {.path {output_dir}/relative_efficiency_{SELECTED_PARAMETER}.png}")

# Print table to console
print(gt_table)
