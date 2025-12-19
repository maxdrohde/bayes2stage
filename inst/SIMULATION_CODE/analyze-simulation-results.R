args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3L) {
    cli::cli_abort(
        "Usage: Rscript analyze-simulation-results.R <parameter_name> <true_value> <sim_setting>"
    )
}

selected_parameter <- args[[1]]
true_param <- as.numeric(args[[2]])
selected_setting <- as.integer(args[[3]])

suppressWarnings(suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
}))

source("simulation_config.R")

theme_set(
    cowplot::theme_cowplot(font_size = 12, font_family = "Source Sans Pro") +
        theme(
            plot.title = element_text(face = "bold"),
            geom = element_geom(pointsize = 3, linewidth = 0.8)
        )
)

# Forest plot with point estimates and CIs
# Expects columns: type, estimate, lower, upper
forest_plot <- function(data, title, x_lab, add_vline = TRUE) {
    result <-
        data |>
        ggplot(aes(y = type, x = estimate)) +
        (if (add_vline) {
            geom_vline(
                xintercept = 0,
                linetype = "dashed",
                color = "gray",
                alpha = 0.7
            )
        }) +
        geom_errorbar(
            aes(xmin = lower, xmax = upper),
            width = 0.2,
            orientation = "y"
        ) +
        geom_point() +
        labs(title = title, x = x_lab, y = "")
    return(result)
}

format_ci <- function(est, lower, upper, digits = 2) {
    result <- dplyr::if_else(
        is.na(est) | is.na(lower) | is.na(upper),
        "NA (NA, NA)",
        paste0(round(est, digits), " (", round(lower, digits), ", ", round(upper, digits), ")")
    )
    return(result)
}

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

cached_file <- "processed_data/combined_data.parquet"

if (!fs::file_exists(cached_file)) {
    cli::cli_abort(
        "Data not found at {.path {cached_file}}. Run batch.R first."
    )
}

full_data <- arrow::read_parquet(cached_file)

# ------------------------------------------------------------------------------
# Filter to selected parameter and sim_setting
# ------------------------------------------------------------------------------

sim_results <- full_data |>
    dplyr::filter(
        parameter == selected_parameter,
        sim_setting == selected_setting
    )

if (nrow(sim_results) == 0) {
    cli::cli_abort(
        "No data found for {selected_parameter} in sim_setting: {selected_setting}"
    )
}

# ------------------------------------------------------------------------------
# MCMC convergence diagnostics
# ------------------------------------------------------------------------------

has_mcmc_diagnostics <- SHOW_MCMC_DIAGNOSTICS &&
    "rhat" %in% names(sim_results) &&
    any(!is.na(sim_results$rhat))

has_hmc_diagnostics <- has_mcmc_diagnostics &&
    "divergent_transitions" %in% names(sim_results) &&
    any(!is.na(sim_results$divergent_transitions))

mcmc_summary <- NULL
hmc_summary <- NULL
p_rhat <- NULL
p_ess <- NULL

if (has_mcmc_diagnostics) {
    mcmc_diagnostics <- sim_results |>
        dplyr::filter(!is.na(rhat)) |>
        dplyr::select(type, rhat, ess_bulk, ess_tail)

    mcmc_summary <- mcmc_diagnostics |>
        dplyr::summarise(
            `Rhat (mean)` = round(mean(rhat, na.rm = TRUE), 3),
            `Rhat (max)` = round(max(rhat, na.rm = TRUE), 3),
            `% Rhat>threshold` = round(
                100 * mean(rhat > DEFAULT_RHAT_THRESHOLD, na.rm = TRUE),
                1
            ),
            `ESS bulk (med)` = round(median(ess_bulk, na.rm = TRUE), 0),
            `ESS bulk (min)` = round(min(ess_bulk, na.rm = TRUE), 0),
            `ESS tail (med)` = round(median(ess_tail, na.rm = TRUE), 0),
            `ESS tail (min)` = round(min(ess_tail, na.rm = TRUE), 0),
            .by = type
        ) |>
        dplyr::rename(Type = type)

    rhat_max <- max(mcmc_diagnostics$rhat, na.rm = TRUE)
    rhat_xlim_upper <- max(1.02, rhat_max * 1.01)

    p_rhat <-
        ggplot(mcmc_diagnostics) +
        aes(x = rhat, y = type) +
        geom_boxplot(
            aes(fill = type),
            alpha = 0.3,
            outlier.shape = NA,
            width = 0.5
        ) +
        geom_jitter(aes(color = type), height = 0.2, alpha = 0.3, size = 0.5) +
        geom_vline(
            xintercept = DEFAULT_RHAT_THRESHOLD,
            linetype = "dashed",
            color = "red",
            linewidth = 0.8
        ) +
        labs(
            title = "Rhat Distribution",
            caption = glue::glue(
                "Dashed line: Rhat = {DEFAULT_RHAT_THRESHOLD} (convergence threshold)"
            ),
            x = "Rhat",
            y = ""
        ) +
        theme(legend.position = "none") +
        coord_cartesian(xlim = c(0.99, rhat_xlim_upper))

    ess_data <- mcmc_diagnostics |>
        tidyr::pivot_longer(
            cols = c(ess_bulk, ess_tail),
            names_to = "ess_type",
            values_to = "ess"
        ) |>
        dplyr::filter(is.finite(ess)) |>
        dplyr::mutate(
            ess_type = dplyr::case_match(
                ess_type,
                "ess_bulk" ~ "Bulk",
                "ess_tail" ~ "Tail"
            )
        )

    p_ess <- ggplot(ess_data, aes(x = ess, y = type, fill = type)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA) +
        geom_jitter(aes(color = type), height = 0.2, alpha = 0.3, size = 0.5) +
        geom_vline(
            xintercept = DEFAULT_ESS_THRESHOLD,
            linetype = "dashed",
            color = "red",
            linewidth = 0.8
        ) +
        facet_wrap(~ess_type) +
        labs(
            title = "Effective Sample Size (ESS)",
            caption = glue::glue(
                "Dashed line: ESS = {DEFAULT_ESS_THRESHOLD} (minimum recommended)"
            ),
            x = "ESS",
            y = ""
        ) +
        theme(legend.position = "none")

    if (has_hmc_diagnostics) {
        hmc_diagnostics <- sim_results |>
            dplyr::filter(!is.na(divergent_transitions))

        hmc_summary <- hmc_diagnostics |>
            dplyr::summarise(
                `Percent with Divergences` =
                    (100 * mean(divergent_transitions > 0, na.rm = TRUE)) |>
                    round(1),
                `mean(Divergences)` =
                    mean(divergent_transitions, na.rm = TRUE) |>
                    round(1),
                `Percent with Max Treedepth` =
                    (100 * mean(max_treedepth_exceeded > 0, na.rm = TRUE)) |>
                    round(1),
                `mean(E-BFMI)` =
                    mean(ebfmi_min, na.rm = TRUE) |>
                    round(2),
                `min(E-BFMI)` =
                    min(ebfmi_min, na.rm = TRUE) |>
                    round(2),
                .by = type
            ) |>
            dplyr::rename(Type = type)
    }
}

# ------------------------------------------------------------------------------
# Paired bootstrap analysis
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Find complete iterations (present in all types for valid pairing)
# ------------------------------------------------------------------------------

all_types <- unique(sim_results$type)
n_types_total <- length(all_types)

iter_type_counts <- sim_results |>
    dplyr::summarise(n_types = dplyr::n_distinct(type), .by = sim_iter) |>
    dplyr::mutate(complete = n_types == n_types_total)

complete_iters <- iter_type_counts$sim_iter[iter_type_counts$complete]

if (!all(iter_type_counts$complete)) {
    incomplete <- dplyr::filter(iter_type_counts, !complete)
    missing_types <- sim_results |>
        dplyr::filter(sim_iter %in% incomplete$sim_iter) |>
        dplyr::summarise(present = list(unique(type)), .by = sim_iter) |>
        dplyr::pull(present) |>
        purrr::map(\(x) setdiff(all_types, x)) |>
        unlist() |>
        table() |>
        sort(decreasing = TRUE)

    cli::cli_alert_warning("Dropped {.val {nrow(incomplete)}} incomplete iterations")
    cli::cli_alert_info("Missing types:")
    cli::cli_bullets(stats::setNames(
        glue::glue("{names(missing_types)}: {missing_types}"),
        rep("*", length(missing_types))
    ))
    cli::cli_alert_info("Using {.val {length(complete_iters)}} complete iterations")
}

if (length(complete_iters) == 0L) {
    cli::cli_abort("No complete iterations found")
}

sim_results_paired <- sim_results |>
    dplyr::filter(sim_iter %in% complete_iters)

n_iters <- length(complete_iters)

# ------------------------------------------------------------------------------
# Paired bootstrap: resample by sim_iter to preserve within-iteration pairing
# Uses data.table for efficient keyed joins and aggregation
# ------------------------------------------------------------------------------

dt <- data.table::as.data.table(sim_results_paired, key = "sim_iter")
iters <- unique(dt$sim_iter)

compute_stats <- function(x, by_cols, truth) {
    x[, .(
        bias     = mean(estimate, na.rm = TRUE) - truth,
        se_emp   = sd(estimate, na.rm = TRUE),
        se_model = mean(se, na.rm = TRUE)
    ), by = by_cols
    ][, se_pct_bias := data.table::fifelse(
        is.finite(se_emp) & se_emp > 0,
        100 * (se_model - se_emp) / se_emp, NA_real_)][]
}

# 1. Original Statistics
orig_stats <- compute_stats(dt, "type", true_param) |> tibble::as_tibble()

# 2. Bootstrap Statistics
boot_map <- data.table::data.table(
    .boot = rep(seq_len(DEFAULT_N_BOOT_REPS), each = length(iters)),
    sim_iter = sample(iters, length(iters) * DEFAULT_N_BOOT_REPS, replace = TRUE)
)

boot_samples <- dt[boot_map, on = "sim_iter", allow.cartesian = TRUE] |>
    compute_stats(c(".boot", "type"), true_param) |>
    tibble::as_tibble()

# ------------------------------------------------------------------------------
# Compute percentile CIs
# ------------------------------------------------------------------------------

alpha <- 1 - DEFAULT_CI_LEVEL

boot_cis <- boot_samples |>
    dplyr::summarise(
        bias_lower = quantile(bias, alpha / 2, na.rm = TRUE),
        bias_upper = quantile(bias, 1 - alpha / 2, na.rm = TRUE),
        se_emp_lower = quantile(se_emp, alpha / 2, na.rm = TRUE),
        se_emp_upper = quantile(se_emp, 1 - alpha / 2, na.rm = TRUE),
        se_model_lower = quantile(se_model, alpha / 2, na.rm = TRUE),
        se_model_upper = quantile(se_model, 1 - alpha / 2, na.rm = TRUE),
        se_pct_bias_lower = quantile(se_pct_bias, alpha / 2, na.rm = TRUE),
        se_pct_bias_upper = quantile(se_pct_bias, 1 - alpha / 2, na.rm = TRUE),
        .by = type
    )

# ------------------------------------------------------------------------------
# Build output structures for downstream plotting code
# ------------------------------------------------------------------------------

stats_with_cis <- dplyr::left_join(orig_stats, boot_cis, by = "type")

se_table <- stats_with_cis |>
    dplyr::mutate(
        diff = se_model - se_emp,
        pct_bias = dplyr::if_else(se_emp > 0, 100 * diff / se_emp, NA_real_)
    )

boot_samples_se <- dplyr::select(boot_samples, type, se_emp, se_model)

extract_ci <- function(data, est_col, lower_col, upper_col) {
    dplyr::transmute(
        data,
        type,
        estimate = .data[[est_col]],
        lower = .data[[lower_col]],
        upper = .data[[upper_col]]
    )
}

boot_bias <- extract_ci(stats_with_cis, "bias", "bias_lower", "bias_upper")
boot_se <- extract_ci(stats_with_cis, "se_emp", "se_emp_lower", "se_emp_upper")
boot_se_pct_bias <- extract_ci(
    stats_with_cis,
    "se_pct_bias",
    "se_pct_bias_lower",
    "se_pct_bias_upper"
)

# ------------------------------------------------------------------------------
# Pairwise differences in empirical SE between methods
# Goal: Test whether one method has significantly different variability than another
# ------------------------------------------------------------------------------

types_vec <- sort(unique(as.character(boot_samples$type)))

if (length(types_vec) >= 2) {
    # Reshape bootstrap samples to wide format
    # Each row = one bootstrap replicate, columns = SE for each method
    # This makes computing pairwise differences easy (just subtract columns)
    boot_se_wide <- boot_samples |>
        dplyr::select(type, .boot, se_emp) |>
        tidyr::pivot_wider(names_from = type, values_from = se_emp)

    # Create named vector for easy lookup of original SE estimates
    orig_se <- stats::setNames(orig_stats$se_emp, orig_stats$type)

    # Generate all unique pairs of methods to compare
    # e.g., for types A, B, C: returns list((A,B), (A,C), (B,C))
    all_pairs <- utils::combn(types_vec, 2, simplify = FALSE)

    # For each pair, compute the difference in SE with bootstrap CI
    pairwise_se_diff <- purrr::map_df(all_pairs, \(pair) {
        type1 <- pair[1]
        type2 <- pair[2]

        # Compute SE difference for each bootstrap replicate
        boot_diffs <- boot_se_wide[[type1]] - boot_se_wide[[type2]]

        # Point estimate: difference in original SE values
        point_estimate <- orig_se[type1] - orig_se[type2]

        # Bootstrap percentile CI for the difference
        ci_lower <- quantile(boot_diffs, alpha / 2, na.rm = TRUE)
        ci_upper <- quantile(boot_diffs, 1 - alpha / 2, na.rm = TRUE)

        # Clean up type names (remove newlines) for display
        label <- paste(gsub("\n", " ", c(type1, type2)), collapse = " - ")

        tibble::tibble(
            comparison = label,
            estimate = point_estimate,
            lower = ci_lower,
            upper = ci_upper
        )
    })
} else {
    # Not enough methods to compare - return empty tibble with correct structure
    pairwise_se_diff <- tibble::tibble(
        comparison = character(),
        estimate = numeric(),
        lower = numeric(),
        upper = numeric()
    )
}

# ------------------------------------------------------------------------------
# Win probability: which method most often has the smallest empirical SE
# ------------------------------------------------------------------------------

boot_winners <- boot_samples |>
    dplyr::filter(is.finite(se_emp)) |>
    dplyr::filter(se_emp == min(se_emp), .by = .boot)

# Handle ties: if multiple types tie for min, give fractional wins
boot_winner_counts <- boot_winners |>
    dplyr::add_count(.boot, name = "n_tied") |>
    dplyr::mutate(win_share = 1 / n_tied) |>
    dplyr::summarise(wins = sum(win_share), .by = type)

# Ensure all types are represented (even with 0 wins)
all_types <- unique(boot_samples$type)
win_probs <- tibble::tibble(type = all_types) |>
    dplyr::left_join(boot_winner_counts, by = "type") |>
    dplyr::mutate(wins = dplyr::coalesce(wins, 0))

# Compute Agresti-Coull CIs for win probability (better than Wald at boundaries)
n_boot <- DEFAULT_N_BOOT_REPS
z <- qnorm((1 + DEFAULT_CI_LEVEL) / 2)
win_probs <- win_probs |>
    dplyr::mutate(
        prop = wins / n_boot,
        win_pct = 100 * prop,
        # Agresti-Coull: add z^2/2 successes and z^2/2 failures
        n_tilde = n_boot + z^2,
        p_tilde = (wins + z^2 / 2) / n_tilde,
        se_tilde = sqrt(p_tilde * (1 - p_tilde) / n_tilde),
        lower = 100 * pmax(0, p_tilde - z * se_tilde),
        upper = 100 * pmin(1, p_tilde + z * se_tilde)
    ) |>
    dplyr::select(-n_tilde, -p_tilde, -se_tilde)

# ------------------------------------------------------------------------------
# Plotting
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# SE comparison plot
# ------------------------------------------------------------------------------

range_vals <- range(
    boot_samples_se$se_emp,
    boot_samples_se$se_model,
    na.rm = TRUE,
    finite = TRUE
)

p_se <-
    ggplot() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    stat_ellipse(
        data = boot_samples_se,
        aes(x = se_emp, y = se_model, fill = type),
        geom = "polygon",
        type = "norm",
        level = DEFAULT_CI_LEVEL,
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
        show.legend = FALSE,
        fontface = "bold",
        max.overlaps = Inf
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

p_se_bias_boot <-
    forest_plot(
        boot_se_pct_bias,
        title = "Percent bias in Model-based SE",
        x_lab = "Percent bias (%)"
    )

p_bias <- boot_bias |>
    ggplot(aes(y = type, x = estimate)) +
    geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "gray",
        alpha = 0.7
    ) +
    geom_errorbar(
        aes(xmin = lower, xmax = upper),
        width = 0.2,
        orientation = "y"
    ) +
    geom_point() +
    labs(title = "Bias", y = "")

if (include_pct_bias) {
    p_bias <- p_bias +
        scale_x_continuous(
            name = "Bias",
            sec.axis = sec_axis(
                ~ . * 100 / abs(true_param),
                name = "Percent bias (%)"
            )
        )
} else {
    p_bias <- p_bias + labs(x = "Bias")
}

p_emp_se <- forest_plot(
    boot_se,
    title = "Empirical SE",
    x_lab = "Empirical SE",
    add_vline = FALSE
)

# ------------------------------------------------------------------------------
# Pairwise SE difference plot
# ------------------------------------------------------------------------------

if (nrow(pairwise_se_diff) > 0) {
    p_pairwise_se <- pairwise_se_diff |>
        ggplot(aes(y = comparison, x = estimate)) +
        geom_vline(
            xintercept = 0,
            linetype = "dashed",
            color = "gray",
            alpha = 0.7
        ) +
        geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2) +
        geom_point() +
        labs(
            title = "Pairwise Differences in Empirical SE",
            x = "Difference in Empirical SE",
            y = ""
        )
} else {
    p_pairwise_se <- ggplot() +
        annotate(
            "text",
            x = 0.5,
            y = 0.5,
            label = "Fewer than 2 types\n(no pairwise comparisons)",
            size = 4
        ) +
        theme_void() +
        labs(title = "Pairwise Differences in Empirical SE")
}

# ------------------------------------------------------------------------------
# Win probability plot
# ------------------------------------------------------------------------------

n_types_for_plot <- length(all_types)
p_win_prob <- win_probs |>
    ggplot(aes(y = type, x = win_pct)) +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.2) +
    geom_point() +
    labs(
        title = "Win Probability (Smallest Empirical SE)",
        x = "Win Probability (%)",
        y = ""
    )

# ------------------------------------------------------------------------------
# Summary table
# ------------------------------------------------------------------------------

pct_bias_col <- if (include_pct_bias) {
    boot_bias |>
        dplyr::mutate(dplyr::across(c(estimate, lower, upper), \(x) {
            x * 100 / abs(true_param)
        })) |>
        dplyr::transmute(
            type,
            `Percent Bias` = format_ci(estimate, lower, upper)
        )
} else {
    tibble::tibble(
        type = unique(sim_results_paired$type),
        `Percent Bias` = "N/A"
    )
}

summary_tables <- list(
    dplyr::count(sim_results_paired, type, name = "N_sims"),
    dplyr::transmute(boot_bias, type, Bias = format_ci(estimate, lower, upper)),
    pct_bias_col,
    dplyr::transmute(
        boot_se,
        type,
        `Empirical SE` = format_ci(estimate, lower, upper)
    ),
    dplyr::transmute(
        se_table,
        type,
        `Model-based SE` = format_ci(se_model, se_model_lower, se_model_upper)
    ),
    dplyr::transmute(
        boot_se_pct_bias,
        type,
        `SE Percent Bias` = format_ci(estimate, lower, upper)
    )
)

combined_table <- purrr::reduce(
    summary_tables,
    dplyr::left_join,
    by = "type"
) |>
    dplyr::rename(Type = type)

# ------------------------------------------------------------------------------
# Get N and sampling_fraction for the current setting
# ------------------------------------------------------------------------------
current_grid_row <- dplyr::slice(get_simulation_grid(), selected_setting)
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

p_table <- wrap_elements(table_grob) + ggtitle("Summary Statistics")

if (has_mcmc_diagnostics) {
    mcmc_table_grob <- gridExtra::tableGrob(
        mcmc_summary,
        rows = NULL,
        theme = table_theme
    )
    p_mcmc_table <- wrap_elements(mcmc_table_grob) + ggtitle("MCMC Diagnostics")
}

if (has_hmc_diagnostics) {
    hmc_table_grob <- gridExtra::tableGrob(
        hmc_summary,
        rows = NULL,
        theme = table_theme
    )
    p_hmc_table <- wrap_elements(hmc_table_grob) + ggtitle("HMC Diagnostics")
}

# ------------------------------------------------------------------------------
# Combined plot
# ------------------------------------------------------------------------------

# Determine X distribution label
x_dist_raw <- current_grid_row[["x_dist"]]
x_dist_label <- if (
    x_dist_raw == "binomial" && identical(current_grid_row[["x_size"]], 1L)
) {
    "Bernoulli"
} else {
    stringr::str_to_title(x_dist_raw)
}

# Build combined plot with appropriate diagnostics panels
plot_title <- glue::glue(
    "Parameter: {format_parameter_name(selected_parameter)} (True value: {true_param})",
    "\nSetting: {selected_setting}, N_total: {current_N}, N_sampled: {current_n_sampled}, X: {x_dist_label}"
)
plot_theme <- theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    plot.caption = element_text(size = 9, hjust = 1)
)

# ------------------------------------------------------------------------------
# Settings summary panel
# ------------------------------------------------------------------------------

settings_col1 <- glue::glue(
    "TARGET
  parameter: {selected_parameter}
  true_value: {true_param}
  sim_setting: {selected_setting}

DATA GENERATION
  N: {current_grid_row$N}
  sampling_frac: {current_grid_row$sampling_fraction}
  M: {current_grid_row$M}
  x_dist: {current_grid_row$x_dist}

SAMPLING DESIGN
  type: {current_grid_row$sampling_type}
  cutoffs: ({DEFAULT_CUTOFF_LOW}, {DEFAULT_CUTOFF_HIGH})
  props: ({DEFAULT_PROP_LOW}, {DEFAULT_PROP_MIDDLE}, {DEFAULT_PROP_HIGH})"
)

settings_col2 <- glue::glue(
    "MODEL COEFFICIENTS
  beta_x: {current_grid_row$beta_x}
  beta_z: {current_grid_row$beta_z}
  beta_t: {current_grid_row$beta_t}
  beta_x:t: {current_grid_row$beta_t_x_interaction}
  beta_z:t: {current_grid_row$beta_t_z_interaction}
  alpha: {current_grid_row$alpha_main}
  sigma: {current_grid_row$error_sd}

RANDOM EFFECTS
  sigma_int: {current_grid_row$rand_intercept_sd}
  sigma_slope: {current_grid_row$rand_slope_sd}
  corr: {current_grid_row$rand_eff_corr}"
)

settings_col3 <- glue::glue(
    "IMPUTATION MODEL
  gamma0: {current_grid_row$gamma0}
  gamma1: {current_grid_row$gamma1}
  gamma2: {current_grid_row$gamma2}
  gamma_sd: {current_grid_row$gamma_sd}

STAN OPTIONS
  distribution: {DEFAULT_STAN_DISTRIBUTION}
  parameterization: {DEFAULT_PARAMETERIZATION}
  inference: {DEFAULT_INFERENCE_ARGS$inference_method}

ANALYSIS
  bootstrap_reps: {DEFAULT_N_BOOT_REPS}
  n_iters: {n_iters}
  generated: {format(Sys.time(), '%Y-%m-%d %H:%M')}"
)

settings_grobs <- purrr::map(
    list(
        settings_col1,
        settings_col2,
        settings_col3
    ),
    \(txt) {
        grid::textGrob(
            txt,
            x = 0.05,
            y = 0.95,
            hjust = 0,
            vjust = 1,
            gp = grid::gpar(fontsize = 9, fontfamily = "mono")
        )
    }
)

p_settings_inner <- wrap_plots(
    wrap_elements(settings_grobs[[1]]),
    wrap_elements(settings_grobs[[2]]),
    wrap_elements(settings_grobs[[3]]),
    nrow = 1
) &
    theme(plot.margin = margin(2, 10, 2, 10))

p_settings <- p_settings_inner +
    plot_annotation(title = "Simulation Settings") &
    theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0),
        plot.background = element_rect(fill = "gray98", color = NA)
    )

if (has_hmc_diagnostics) {
    combined_plot <-
        wrap_plots(
            p_table,
            p_bias,
            p_se,
            p_se_bias_boot,
            p_emp_se,
            p_pairwise_se,
            p_win_prob,
            p_rhat,
            p_ess,
            p_mcmc_table,
            p_hmc_table,
            p_settings,
            ncol = 1,
            heights = c(1, 1.5, 2, 1.5, 1.5, 2, 1.5, 1, 1, 1.5, 1.5, 1.5)
        ) +
        plot_annotation(title = plot_title, theme = plot_theme)
} else if (has_mcmc_diagnostics) {
    combined_plot <-
        wrap_plots(
            p_table,
            p_bias,
            p_se,
            p_se_bias_boot,
            p_emp_se,
            p_pairwise_se,
            p_win_prob,
            p_rhat,
            p_ess,
            p_mcmc_table,
            p_settings,
            ncol = 1,
            heights = c(1, 1.5, 2, 1.5, 1.5, 2, 1.5, 1, 1, 1.5, 1.5)
        ) +
        plot_annotation(title = plot_title, theme = plot_theme)
} else {
    combined_plot <-
        wrap_plots(
            p_table,
            p_bias,
            p_se,
            p_se_bias_boot,
            p_emp_se,
            p_pairwise_se,
            p_win_prob,
            p_settings,
            ncol = 1,
            heights = c(2, 0.8, 2, 0.8, 0.8, 0.8, 0.8, 1.5)
        ) +
        plot_annotation(title = plot_title, theme = plot_theme)
}

# ------------------------------------------------------------------------------
# Save output
# ------------------------------------------------------------------------------
plots_dir <- glue::glue("plots/setting_{selected_setting}")
fs::dir_create(plots_dir)

ggsave(
    plot = combined_plot,
    filename = glue::glue(
        "{plots_dir}/simulation_plots_{selected_parameter}.png"
    ),
    units = "in",
    width = 14,
    height = 20,
    device = ragg::agg_png,
    bg = "white",
    dpi = 300
)
