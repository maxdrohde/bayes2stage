# ==============================================================================
# CONFIGURATION
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2L) {
  stop("Usage: Rscript analyze-simulation-results.R <parameter_name> <true_value> [sim_setting]")
}

selected_parameter <- args[[1]]
true_param <- as.numeric(args[[2]])

# Optional sim_setting argument (3rd arg) - defaults to NULL for backward compat
selected_setting <- if (length(args) >= 3L) as.integer(args[[3]]) else NULL

# ==============================================================================
# LOAD LIBRARIES
# ==============================================================================

suppressPackageStartupMessages(library(tidyverse))

library(patchwork)

source("simulation_config.R")

theme_set(
  theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      geom = element_geom(
        pointsize = 3,
        linewidth = 0.8
      )
    )
)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# clean_type_names(), format_parameter_name(), and TYPE_LEVELS
# are defined in simulation_config.R

safe_boot_ci <- function(boot_obj, index = 1) {
  if (is.na(index)) {
    return(c(NA_real_, NA_real_))
  }

  if (any(!is.finite(boot_obj$t[, index]))) {
    return(c(NA_real_, NA_real_))
  }
  tryCatch(
    {
      ci <- boot::boot.ci(boot_obj, type = "bca", index = index)$bca
      c(ci[4], ci[5])
    },
    error = function(e) c(NA_real_, NA_real_)
  )
}

plot_bootstrap_metrics <- function(data,
                                   title,
                                   x_lab,
                                   add_vline = TRUE) {
  p <-
    data |>
    ggplot(aes(
      y = type,
      x = estimate
    ))

  if (add_vline) {
    p <- p + geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "gray50",
      alpha = 0.7
    )
  }

  p +
    geom_errorbar(
      aes(
        xmin = lower,
        xmax = upper
      ),
      width = 0.2,
      orientation = "y"
    ) +
    geom_point() +
    labs(
      title = title,
      x = x_lab,
      y = ""
    )
}

format_ci <- function(est,
                      lower,
                      upper,
                      digits = 2) {
  if (any(is.na(c(est, lower, upper)))) {
    return("NA (NA, NA)")
  }
  glue::glue("{round(est, digits)} ({round(lower, digits)}, {round(upper, digits)})")
}

# ==============================================================================
# LOAD DATA
# ==============================================================================

# Load cached data created by batch.R
cached_file <- "processed_data/combined_data.parquet"

if (!fs::file_exists(cached_file)) {
  stop(
    "Cached data not found at ", cached_file,
    ". Run batch.R first to merge parquet files."
  )
}

full_data <- arrow::read_parquet(cached_file)

# ------------------------------------------------------------------------------
# Filter to selected parameter (and optionally sim_setting)
# ------------------------------------------------------------------------------

res_x <-
  full_data |>
  filter(parameter == selected_parameter)

# Filter by sim_setting if provided
if (!is.null(selected_setting)) {
  # Handle backward compatibility: if sim_setting column doesn't exist, treat all as setting 1
  if ("sim_setting" %in% names(res_x)) {
    res_x <- res_x |> filter(sim_setting == selected_setting)
  } else if (selected_setting != 1L) {
    # No sim_setting column and not requesting setting 1 - no data available
    res_x <- res_x[0, ]
  }
  # If selected_setting == 1 and no column exists, keep all data (backward compat)
}

if (nrow(res_x) == 0) {
  setting_msg <- if (!is.null(selected_setting)) glue::glue(", sim_setting: {selected_setting}") else ""
  stop(glue::glue("No data found for parameter: {selected_parameter}{setting_msg}"))
}

# ==============================================================================
# MCMC CONVERGENCE DIAGNOSTICS
# ==============================================================================

has_mcmc_diagnostics <- SHOW_MCMC_DIAGNOSTICS &&
    "rhat" %in% names(res_x) &&
    any(!is.na(res_x$rhat))

has_hmc_diagnostics <- has_mcmc_diagnostics &&
    "divergent_transitions" %in% names(res_x) &&
    any(!is.na(res_x$divergent_transitions))

if (has_mcmc_diagnostics) {
    mcmc_summary <- res_x |>
        filter(!is.na(rhat)) |>
        group_by(type) |>
        summarise(
            `Rhat (mean)` = round(mean(rhat, na.rm = TRUE), 3),
            `Rhat (max)` = round(max(rhat, na.rm = TRUE), 3),
            `% Rhat>1.01` = round(100 * mean(rhat > 1.01, na.rm = TRUE), 1),
            `ESS bulk (med)` = round(median(ess_bulk, na.rm = TRUE), 0),
            `ESS bulk (min)` = round(min(ess_bulk, na.rm = TRUE), 0),
            `ESS tail (med)` = round(median(ess_tail, na.rm = TRUE), 0),
            `ESS tail (min)` = round(min(ess_tail, na.rm = TRUE), 0),
            .groups = "drop"
        ) |>
        rename(Type = type)

    # Create MCMC diagnostic plots
    mcmc_data <- res_x |>
        filter(!is.na(rhat)) |>
        select(type, rhat, ess_bulk, ess_tail)

    p_rhat <- ggplot(mcmc_data, aes(x = rhat, y = type, fill = type)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        geom_jitter(aes(color = type), height = 0.2, alpha = 0.3, size = 0.5) +
        geom_vline(xintercept = 1.01, linetype = "dashed", color = "red", linewidth = 0.8) +
        labs(title = "Rhat Distribution", x = "Rhat", y = "") +
        theme(legend.position = "none") +
        coord_cartesian(xlim = c(0.99, max(1.02, max(mcmc_data$rhat, na.rm = TRUE) * 1.01)))

    p_ess <- mcmc_data |>
        tidyr::pivot_longer(cols = c(ess_bulk, ess_tail), names_to = "ess_type", values_to = "ess") |>
        filter(is.finite(ess)) |>
        mutate(ess_type = ifelse(ess_type == "ess_bulk", "Bulk", "Tail")) |>
        ggplot(aes(x = type, y = ess, fill = type)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        geom_jitter(aes(color = type), width = 0.2, alpha = 0.3, size = 0.5) +
        geom_hline(yintercept = 400, linetype = "dashed", color = "orange", linewidth = 0.8) +
        facet_wrap(~ess_type) +
        labs(title = "Effective Sample Size", x = "", y = "ESS") +
        theme(legend.position = "none")

    # HMC-specific diagnostics (divergences, treedepth, E-BFMI)
    if (has_hmc_diagnostics) {
        hmc_summary <- res_x |>
            filter(!is.na(divergent_transitions)) |>
            group_by(type) |>
            summarise(
                `% w/ Divergences` = round(100 * mean(divergent_transitions > 0, na.rm = TRUE), 1),
                `Divergences (mean)` = round(mean(divergent_transitions, na.rm = TRUE), 1),
                `% Max Treedepth` = round(100 * mean(max_treedepth_exceeded > 0, na.rm = TRUE), 1),
                `E-BFMI (mean)` = round(mean(ebfmi_min, na.rm = TRUE), 2),
                `E-BFMI (min)` = round(min(ebfmi_min, na.rm = TRUE), 2),
                .groups = "drop"
            ) |>
            rename(Type = type)

        # Bar chart of HMC issues
        hmc_issue_data <- res_x |>
            filter(!is.na(divergent_transitions)) |>
            group_by(type) |>
            summarise(
                `Divergent` = mean(divergent_transitions > 0, na.rm = TRUE),
                `Max Treedepth` = mean(max_treedepth_exceeded > 0, na.rm = TRUE),
                `Low E-BFMI` = mean(ebfmi_min < 0.3, na.rm = TRUE),
                .groups = "drop"
            ) |>
            tidyr::pivot_longer(cols = -type, names_to = "issue", values_to = "proportion")

        p_hmc <- ggplot(hmc_issue_data, aes(x = type, y = proportion * 100, fill = issue)) +
            geom_col(position = "dodge", alpha = 0.8) +
            geom_hline(yintercept = 0, color = "gray50") +
            labs(
                title = "HMC Diagnostic Issues",
                x = "",
                y = "% of Datasets",
                fill = "Issue Type"
            ) +
            scale_fill_brewer(palette = "Set2") +
            theme(legend.position = "bottom")
    } else {
        hmc_summary <- NULL
        p_hmc <- NULL
    }
} else {
    mcmc_summary <- NULL
    hmc_summary <- NULL
    p_rhat <- NULL
    p_ess <- NULL
    p_hmc <- NULL
}

# ==============================================================================
# BOOTSTRAP ANALYSIS
# ==============================================================================

# ------------------------------------------------------------------------------
# SE comparison (empirical vs model-based)
# ------------------------------------------------------------------------------

combined_stat_fun <- function(data, indices) {
  d <- data[indices, ]
  se_emp <- sd(d$estimate, na.rm = TRUE)
  se_model <- mean(d$se, na.rm = TRUE)
  bias <- mean(d$estimate, na.rm = TRUE) - true_param
  se_pct_bias <- 100 * (se_model - se_emp) / se_emp

  c(
    bias = bias,
    se_emp = se_emp,
    se_model = se_model,
    se_pct_bias = se_pct_bias
  )
}

combined_boots <- res_x |>
  group_split(type, .keep = TRUE) |>
  map_df(function(d) {
    type_name <- unique(d$type)
    boot_obj <- boot::boot(
      data = select(d, -type),
      statistic = combined_stat_fun,
      R = DEFAULT_N_BOOT_REPS
    )

    stat_names <- names(boot_obj$t0)
    if (is.null(stat_names) || length(stat_names) == 0) {
      stat_names <- paste0("stat", seq_along(boot_obj$t0))
    }
    colnames(boot_obj$t) <- stat_names

    tibble(
      type = type_name,
      boot_obj = list(boot_obj)
    )
  })

se_table <- combined_boots |>
  mutate(
    se_emp = map_dbl(boot_obj, ~ .x$t0[["se_emp"]]),
    se_model = map_dbl(boot_obj, ~ .x$t0[["se_model"]]),
    se_model_ci = map(boot_obj, ~ safe_boot_ci(.x, index = match("se_model", names(.x$t0)))),
    se_model_lower = map_dbl(se_model_ci, 1),
    se_model_upper = map_dbl(se_model_ci, 2),
    diff = se_model - se_emp,
    pct_bias = 100 * diff / se_emp
  ) |>
  select(-boot_obj, -se_model_ci)

boot_samples_se <- combined_boots |>
  mutate(
    samples = map(boot_obj, ~ tibble(se_emp = .x$t[, "se_emp"], se_model = .x$t[, "se_model"]))
  ) |>
  select(type, samples) |>
  unnest(samples)

# ==============================================================================
# PLOTTING
# ==============================================================================

# ------------------------------------------------------------------------------
# SE comparison plot
# ------------------------------------------------------------------------------

range_vals <- range(c(boot_samples_se$se_emp, boot_samples_se$se_model))

p_se <-
  ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  stat_ellipse(
    data = boot_samples_se,
    aes(x = se_emp, y = se_model, fill = type),
    geom = "polygon",
    type = "norm",
    level = 0.95,
    color = NA,
    alpha = 0.2
  ) +
  geom_point(
    data = se_table,
    aes(x = se_emp, y = se_model, color = type)
  ) +
  ggrepel::geom_text_repel(
    data = se_table,
    aes(x = se_emp, y = se_model, label = type, color = type),
    show.legend = FALSE, fontface = "bold", max.overlaps = Inf
  ) +
  coord_fixed(xlim = range_vals, ylim = range_vals, clip = "off") +
  labs(
    title = "Empirical SE vs Model-based SE",
    x = "Empirical SE",
    y = "Model-based SE"
  ) +
  theme(
    legend.position = "none",
    plot.margin = margin(10, 10, 20, 10) # add space for labels drawn outside the panel
  )

include_pct_bias <- !isTRUE(all.equal(true_param, 0))

extract_metric <- function(df, metric) {
  df |>
    transmute(
      type,
      estimate = map_dbl(boot_obj, ~ .x$t0[[metric]]),
      ci = map(boot_obj, ~ safe_boot_ci(.x, index = match(metric, names(.x$t0)))),
      lower = map_dbl(ci, 1),
      upper = map_dbl(ci, 2)
    )
}

boot_bias <- extract_metric(combined_boots, "bias")
boot_se <- extract_metric(combined_boots, "se_emp")
boot_results <- extract_metric(combined_boots, "se_pct_bias")

p_se_bias_boot <-
  plot_bootstrap_metrics(boot_results,
    title = "Percent bias in Model-based SE",
    x_lab = "Percent bias (%)"
  )

p_1 <- boot_bias |>
  ggplot(aes(y = type, x = estimate)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2, orientation = "y") +
  geom_point() +
  labs(title = "Bias", y = "")

if (include_pct_bias) {
  p_1 <- p_1 + scale_x_continuous(
    name = "Bias",
    sec.axis = sec_axis(~ . * 100 / abs(true_param), name = "Percent bias (%)")
  )
} else {
  p_1 <- p_1 + labs(x = "Bias")
}

p_3 <-
  plot_bootstrap_metrics(boot_se,
    title = "Empirical SE",
    x_lab = "Empirical SE",
    add_vline = FALSE
  )

# ==============================================================================
# SUMMARY TABLE
# ==============================================================================

count_table <- res_x |>
  count(type, name = "n")

bias_table <- boot_bias |>
  transmute(type, Bias = format_ci(estimate, lower, upper))

pct_bias_table <-
  if (include_pct_bias) {
    boot_bias |>
      mutate(
        pct_est = estimate * 100 / abs(true_param),
        pct_lower = lower * 100 / abs(true_param),
        pct_upper = upper * 100 / abs(true_param)
      ) |>
      transmute(type, `Percent Bias` = format_ci(pct_est, pct_lower, pct_upper))
  } else {
    tibble(
      type = unique(res_x$type),
      `Percent Bias` = "N/A"
    )
  }

emp_se_table <- boot_se |>
  transmute(type, `Empirical SE` = format_ci(estimate, lower, upper))

model_se_table <- se_table |>
  transmute(type, `Model-based SE` = format_ci(se_model, se_model_lower, se_model_upper))

se_pct_bias_table <- boot_results |>
  transmute(type, `SE Percent Bias` = format_ci(estimate, lower, upper))

combined_table <- count_table |>
  left_join(bias_table, by = "type") |>
  left_join(pct_bias_table, by = "type") |>
  left_join(emp_se_table, by = "type") |>
  left_join(model_se_table, by = "type") |>
  left_join(se_pct_bias_table, by = "type") |>
  mutate(Type = type, N_sims = n) |>
  select(Type, N_sims, Bias, `Percent Bias`, `Empirical SE`, `Model-based SE`, `SE Percent Bias`)

# ------------------------------------------------------------------------------
# Get N and sampling_fraction for the current setting
# ------------------------------------------------------------------------------
current_grid_row <- dplyr::slice(get_simulation_grid(), selected_setting %||% 1)
current_N <- current_grid_row[["N"]]
current_sampled_frac <- current_grid_row[["sampling_fraction"]]
current_n_sampled <- as.integer(current_N * current_sampled_frac)

table_theme <- gridExtra::ttheme_default(
  base_size = 11,
  base_family = "Helvetica",
  padding = grid::unit(c(4, 3), "mm"),
  core = list(
    fg_params = list(hjust = 0.5, x = 0.5),
    bg_params = list(fill = c("white", "gray95"))
  ),
  colhead = list(
    fg_params = list(fontface = "bold", hjust = 0.5, x = 0.5),
    bg_params = list(fill = "gray90", col = "white")
  )
)

table_grob <- gridExtra::tableGrob(
  combined_table,
  rows = NULL,
  theme = table_theme
)

p_table <- wrap_elements(table_grob)

# Create MCMC diagnostics table (if available)
if (has_mcmc_diagnostics) {
    mcmc_table_grob <- gridExtra::tableGrob(
        mcmc_summary,
        rows = NULL,
        theme = table_theme
    )
    p_mcmc_table <- wrap_elements(mcmc_table_grob)
}

# Create HMC diagnostics table (if available)
if (has_hmc_diagnostics) {
    hmc_table_grob <- gridExtra::tableGrob(
        hmc_summary,
        rows = NULL,
        theme = table_theme
    )
    p_hmc_table <- wrap_elements(hmc_table_grob)
}

# ==============================================================================
# COMBINED PLOT
# ==============================================================================

# Determine X distribution label
sim_grid <- get_simulation_grid()
x_dist_raw <- sim_grid$x_dist[1]
x_dist_label <- if (x_dist_raw == "binomial" && sim_grid$x_size[1] == 1) {
  "Bernoulli"
} else {
  stringr::str_to_title(x_dist_raw)
}

# Build combined plot with appropriate diagnostics panels
plot_title <- glue::glue(
    "Parameter: {format_parameter_name(selected_parameter)} (True value: {true_param})",
    "{if (!is.null(selected_setting)) paste0('\nSetting: ', selected_setting) else ''}",
    "\nN_total: {current_N}, N_sampled: {current_n_sampled}, X: {x_dist_label}"
)
plot_theme <- theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.caption = element_text(size = 9, hjust = 1)
)

if (has_hmc_diagnostics) {
    # Full layout with MCMC + HMC diagnostics
    combined_plot <-
        wrap_plots(
            A = p_1,
            B = p_3,
            C = p_se,
            D = p_se_bias_boot,
            E = p_table,
            F = p_rhat,
            G = p_ess,
            H = p_mcmc_table,
            I = p_hmc,
            J = p_hmc_table,
            design = "
AABBB
CCDDD
CCEEE
FFGGG
FFHHH
IIIII
JJJJJ
"
        ) +
        plot_annotation(title = plot_title, theme = plot_theme)
} else if (has_mcmc_diagnostics) {
    # MCMC diagnostics only (no HMC)
    combined_plot <-
        wrap_plots(
            A = p_1,
            B = p_3,
            C = p_se,
            D = p_se_bias_boot,
            E = p_table,
            F = p_rhat,
            G = p_ess,
            H = p_mcmc_table,
            design = "
AABBB
CCDDD
CCEEE
FFGGG
FFHHH
"
        ) +
        plot_annotation(title = plot_title, theme = plot_theme)
} else {
    # No MCMC diagnostics
    combined_plot <-
        wrap_plots(
            A = p_1,
            B = p_3,
            C = p_se,
            D = p_se_bias_boot,
            E = p_table,
            design = "
AABBB
CCDDD
CCEEE
"
        ) +
        plot_annotation(title = plot_title, theme = plot_theme)
}

# ==============================================================================
# SAVE OUTPUT
# ==============================================================================
plots_dir <- if (!is.null(selected_setting)) {
  glue::glue("plots/setting_{selected_setting}")
} else {
  "plots"
}
fs::dir_create(plots_dir)

plot_height <- if (has_hmc_diagnostics) {
    19  # Full MCMC + HMC diagnostics
} else if (has_mcmc_diagnostics) {
    15  # MCMC diagnostics only
} else {
    9   # No diagnostics
}

ggsave(
  plot = combined_plot,
  filename = glue::glue("{plots_dir}/simulation_plots_{selected_parameter}.png"),
  units = "in",
  height = plot_height,
  width = 16,
  device = ragg::agg_png
)

# ==============================================================================
# CLEANUP
# ==============================================================================
