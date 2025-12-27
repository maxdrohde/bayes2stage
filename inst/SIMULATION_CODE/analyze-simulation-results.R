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
    cowplot::theme_cowplot(font_size = 11, font_family = "Source Sans Pro") +
        theme(
            plot.title = element_text(face = "bold"),
            geom = element_geom(pointsize = 3, linewidth = 0.8)
        )
)

# Forest plot with point estimates and CIs
# Expects columns: type, estimate, lower, upper
forest_plot <- function(data, title, x_lab, add_vline = FALSE) {
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
        paste0(
            round(est, digits),
            " (",
            round(lower, digits),
            ", ",
            round(upper, digits),
            ")"
        )
    )
    return(result)
}

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------

cached_dataset <- "processed_data/combined_data"

if (!fs::dir_exists(cached_dataset)) {
    cli::cli_abort(
        "Data not found at {.path {cached_dataset}}. Run batch.R first."
    )
}

# Use Arrow predicate pushdown on partitioned dataset: reads only the relevant partition
# Each parallel worker only loads its specific (parameter, sim_setting) subset
sim_results <- arrow::open_dataset(cached_dataset) |>
    dplyr::filter(
        parameter == selected_parameter,
        sim_setting == selected_setting
    ) |>
    dplyr::collect()

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
    rhat_xlim_upper <- max(DEFAULT_RHAT_THRESHOLD + 0.01, rhat_max * 1.01)

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
                `Percent with Divergences` = (100 *
                    mean(divergent_transitions > 0, na.rm = TRUE)) |>
                    round(1),
                `mean(Divergences)` = mean(
                    divergent_transitions,
                    na.rm = TRUE
                ) |>
                    round(1),
                `Percent with Max Treedepth` = (100 *
                    mean(max_treedepth_exceeded > 0, na.rm = TRUE)) |>
                    round(1),
                `mean(E-BFMI)` = mean(ebfmi_min, na.rm = TRUE) |>
                    round(2),
                `min(E-BFMI)` = min(ebfmi_min, na.rm = TRUE) |>
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

    cli::cli_alert_warning(
        "Dropped {.val {nrow(incomplete)}} incomplete iterations"
    )
    cli::cli_alert_info("Missing types:")
    cli::cli_bullets(stats::setNames(
        glue::glue("{names(missing_types)}: {missing_types}"),
        rep("*", length(missing_types))
    ))
    cli::cli_alert_info(
        "Using {.val {length(complete_iters)}} complete iterations"
    )
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
    x[,
        .(
            bias = mean(estimate, na.rm = TRUE) - truth,
            se_emp = sd(estimate, na.rm = TRUE),
            se_model = mean(se, na.rm = TRUE)
        ),
        by = by_cols
    ][,
        se_pct_bias := data.table::fifelse(
            is.finite(se_emp) & se_emp > 0,
            100 * (se_model - se_emp) / se_emp,
            NA_real_
        )
    ][]
}

compute_accuracy_stats <- function(x, by_cols, truth) {
    x[,
        .(
            rmse = sqrt(mean((estimate - truth)^2, na.rm = TRUE)),
            mae = mean(abs(estimate - truth), na.rm = TRUE),
            median_ae = median(abs(estimate - truth), na.rm = TRUE)
        ),
        by = by_cols
    ]
}

compute_coverage_stats <- function(x, by_cols, truth) {
    x[,
        .(
            coverage = mean(
                truth >= q2_5 & truth <= q97_5,
                na.rm = TRUE
            ),
            ci_width_mean = mean(q97_5 - q2_5, na.rm = TRUE),
            n_covered = sum(
                truth >= q2_5 & truth <= q97_5,
                na.rm = TRUE
            )
        ),
        by = by_cols
    ]
}

# 1. Original Statistics
orig_stats <- compute_stats(dt, "type", true_param) |> tibble::as_tibble()
orig_accuracy <- compute_accuracy_stats(dt, "type", true_param) |> tibble::as_tibble()
orig_coverage <- compute_coverage_stats(dt, "type", true_param) |> tibble::as_tibble()

# 2. Chunked Bootstrap (memory-efficient: processes in batches instead of all at once)
chunk_size <- 500L
n_chunks <- ceiling(DEFAULT_N_BOOT_REPS / chunk_size)

boot_results <- purrr::map(seq_len(n_chunks), \(chunk_idx) {
    start_rep <- (chunk_idx - 1L) * chunk_size + 1L
    end_rep <- min(chunk_idx * chunk_size, DEFAULT_N_BOOT_REPS)
    n_reps <- end_rep - start_rep + 1L

    chunk_map <- data.table::data.table(
        .boot = rep(seq(start_rep, end_rep), each = length(iters)),
        sim_iter = sample(iters, length(iters) * n_reps, replace = TRUE)
    )

    # Single join per chunk, compute all stats
    joined <- dt[chunk_map, on = "sim_iter", allow.cartesian = TRUE]

    list(
        stats = compute_stats(joined, c(".boot", "type"), true_param),
        accuracy = compute_accuracy_stats(joined, c(".boot", "type"), true_param),
        coverage = compute_coverage_stats(joined, c(".boot", "type"), true_param)
    )
}, .progress = TRUE)

# Combine chunks
boot_samples <- purrr::map(boot_results, "stats") |>
    data.table::rbindlist() |>
    tibble::as_tibble()

boot_accuracy <- purrr::map(boot_results, "accuracy") |>
    data.table::rbindlist() |>
    tibble::as_tibble()

boot_coverage <- purrr::map(boot_results, "coverage") |>
    data.table::rbindlist() |>
    tibble::as_tibble() |>
    dplyr::mutate(coverage_pct = 100 * coverage)

rm(boot_results)
invisible(gc())

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

# RMSE & MAE CIs
boot_accuracy_cis <- boot_accuracy |>
    dplyr::summarise(
        rmse_lower = quantile(rmse, alpha / 2, na.rm = TRUE),
        rmse_upper = quantile(rmse, 1 - alpha / 2, na.rm = TRUE),
        mae_lower = quantile(mae, alpha / 2, na.rm = TRUE),
        mae_upper = quantile(mae, 1 - alpha / 2, na.rm = TRUE),
        .by = type
    )

accuracy_with_cis <- dplyr::left_join(
    orig_accuracy,
    boot_accuracy_cis,
    by = "type"
)

# Coverage CIs
boot_coverage_cis <- boot_coverage |>
    dplyr::summarise(
        coverage_lower = quantile(coverage, alpha / 2, na.rm = TRUE),
        coverage_upper = quantile(coverage, 1 - alpha / 2, na.rm = TRUE),
        ci_width_lower = quantile(ci_width_mean, alpha / 2, na.rm = TRUE),
        ci_width_upper = quantile(ci_width_mean, 1 - alpha / 2, na.rm = TRUE),
        .by = type
    )

coverage_with_cis <- dplyr::left_join(
    orig_coverage,
    boot_coverage_cis,
    by = "type"
) |>
    dplyr::mutate(
        coverage_pct = 100 * coverage,
        coverage_lower_pct = 100 * coverage_lower,
        coverage_upper_pct = 100 * coverage_upper,
        coverage_deviation = coverage_pct - (100 * DEFAULT_CI_LEVEL)
    )

# Relative Efficiency with CIs
boot_var_wide <- boot_samples |>
    dplyr::mutate(variance = se_emp^2) |>
    dplyr::select(type, .boot, variance) |>
    tidyr::pivot_wider(names_from = type, values_from = variance)

baseline_type <- DEFAULT_BASELINE_TYPE
orig_var <- stats::setNames((orig_stats$se_emp)^2, orig_stats$type)

rel_eff_results <- purrr::map_df(
    setdiff(names(orig_var), baseline_type),
    \(comp_type) {
        boot_ratios <- boot_var_wide[[baseline_type]] /
            boot_var_wide[[comp_type]]
        point_est <- orig_var[baseline_type] / orig_var[comp_type]

        tibble::tibble(
            type = comp_type,
            rel_eff = point_est,
            rel_eff_lower = quantile(boot_ratios, alpha / 2, na.rm = TRUE),
            rel_eff_upper = quantile(boot_ratios, 1 - alpha / 2, na.rm = TRUE)
        )
    }
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
        plot.margin = margin(5, 5, 10, 5) # add space for labels drawn outside the panel
    )

include_pct_bias <- !isTRUE(all.equal(true_param, 0))

p_se_bias_boot <-
    forest_plot(
        boot_se_pct_bias,
        title = "Percent bias in Model-based SE",
        x_lab = "Percent bias (%)",
        add_vline = TRUE
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
    scale_x_continuous(limits = c(0, 100)) +
    labs(
        title = "Win Probability (Smallest Empirical SE)",
        x = "Win Probability (%)",
        y = ""
    )

# ------------------------------------------------------------------------------
# RMSE Forest Plot
# ------------------------------------------------------------------------------

rmse_forest_data <- accuracy_with_cis |>
    dplyr::transmute(
        type,
        estimate = rmse,
        lower = rmse_lower,
        upper = rmse_upper
    )

p_rmse <- forest_plot(
    rmse_forest_data,
    title = "Root Mean Squared Error (RMSE)",
    x_lab = "RMSE",
    add_vline = FALSE
)

# ------------------------------------------------------------------------------
# MAE Forest Plot
# ------------------------------------------------------------------------------

mae_forest_data <- accuracy_with_cis |>
    dplyr::transmute(
        type,
        estimate = mae,
        lower = mae_lower,
        upper = mae_upper
    )

p_mae <- forest_plot(
    mae_forest_data,
    title = "Mean Absolute Error (MAE)",
    x_lab = "MAE",
    add_vline = FALSE
)

# ------------------------------------------------------------------------------
# Coverage Forest Plot
# ------------------------------------------------------------------------------

p_coverage <- coverage_with_cis |>
    ggplot(aes(y = type, x = coverage_pct)) +
    geom_vline(
        xintercept = 100 * DEFAULT_CI_LEVEL,
        linetype = "solid",
        color = "red",
        linewidth = 0.8
    ) +
    geom_errorbar(
        aes(xmin = coverage_lower_pct, xmax = coverage_upper_pct),
        width = 0.2
    ) +
    geom_point(size = 3) +
    labs(
        title = glue::glue("Coverage of {100 * DEFAULT_CI_LEVEL}% Intervals"),
        x = "Coverage (%)",
        y = "",
        caption = "Red line: nominal coverage\nBayesian: credible intervals | ACML: confidence intervals"
    )

# ------------------------------------------------------------------------------
# CI Width Forest Plot
# ------------------------------------------------------------------------------

ci_width_forest_data <- coverage_with_cis |>
    dplyr::transmute(
        type,
        estimate = ci_width_mean,
        lower = ci_width_lower,
        upper = ci_width_upper
    )

p_ci_width <- forest_plot(
    ci_width_forest_data,
    title = "Mean CI Width",
    x_lab = "CI Width",
    add_vline = FALSE
)

# ------------------------------------------------------------------------------
# Relative Efficiency Forest Plot
# ------------------------------------------------------------------------------

p_rel_eff <- rel_eff_results |>
    ggplot(aes(y = type, x = rel_eff)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
    geom_errorbar(
        aes(xmin = rel_eff_lower, xmax = rel_eff_upper),
        width = 0.2
    ) +
    geom_point(size = 3) +
    scale_x_log10() +
    annotation_logticks(sides = "b") +
    labs(
        title = glue::glue("Relative Efficiency vs {baseline_type}"),
        x = "Relative Efficiency (log scale)",
        y = "",
        caption = ">1: more efficient than baseline | <1: less efficient"
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
            `Percent\nBias` = format_ci(estimate, lower, upper)
        )
} else {
    tibble::tibble(
        type = unique(sim_results_paired$type),
        `Percent\nBias` = "N/A"
    )
}

rmse_table_col <- accuracy_with_cis |>
    dplyr::transmute(
        type,
        `Root Mean\nSq. Error` = format_ci(rmse, rmse_lower, rmse_upper)
    )

mae_table_col <- accuracy_with_cis |>
    dplyr::transmute(
        type,
        `Mean\nAbs. Error` = format_ci(mae, mae_lower, mae_upper)
    )

coverage_table_col <- coverage_with_cis |>
    dplyr::transmute(
        type,
        `Coverage\n(%)` = format_ci(
            coverage_pct,
            coverage_lower_pct,
            coverage_upper_pct,
            digits = 1
        )
    )

summary_tables <- list(
    dplyr::count(sim_results_paired, type, name = "N\nSims"),
    dplyr::transmute(boot_bias, type, Bias = format_ci(estimate, lower, upper)),
    pct_bias_col,
    rmse_table_col,
    mae_table_col,
    coverage_table_col,
    dplyr::transmute(
        boot_se,
        type,
        `Empirical\nSE` = format_ci(estimate, lower, upper)
    ),
    dplyr::transmute(
        se_table,
        type,
        `Model-based\nSE` = format_ci(se_model, se_model_lower, se_model_upper)
    ),
    dplyr::transmute(
        boot_se_pct_bias,
        type,
        `SE Percent\nBias (%)` = format_ci(estimate, lower, upper)
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
    base_size = 10,
    base_family = "Helvetica",
    padding = grid::unit(c(4, 4), "mm"),
    core = list(
        fg_params = list(hjust = 0.5, x = 0.5),
        bg_params = list(
            fill = c("white", "#ECF0F1"),
            col = "white",
            lwd = 1.5
        )
    ),
    colhead = list(
        fg_params = list(
            fontsize = 11,
            fontface = "bold",
            col = "white",
            hjust = 0.5,
            x = 0.5
        ),
        bg_params = list(
            fill = "#2C3E50",
            col = "white",
            lwd = 1.5
        )
    )
)

table_grob <- gridExtra::tableGrob(
    combined_table,
    rows = NULL,
    theme = table_theme
)

p_table <- wrap_elements(table_grob)

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
    "Parameter: {format_parameter_name(selected_parameter)}"
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

# Combine all settings text into one vertical block
settings_text_combined <- paste(
    settings_col1,
    settings_col2,
    settings_col3,
    sep = "\n\n"
)

settings_grob <- grid::textGrob(
    settings_text_combined,
    x = 0.05,
    y = 0.95,
    hjust = 0,
    vjust = 1,
    gp = grid::gpar(fontsize = 11, fontfamily = "Source Code Pro")
)

p_settings <- wrap_elements(settings_grob) +
    plot_annotation(title = "Simulation Settings") &
    theme(
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "gray98", color = NA)
    )


# Organize plots into a unified list
# Plots suitable for the grid
grid_plots <- list(
    p_bias,
    p_rmse,
    p_mae,
    p_coverage,
    p_ci_width,
    p_rel_eff,
    p_se,
    p_se_bias_boot,
    p_emp_se,
    p_pairwise_se,
    p_win_prob
)

if (has_mcmc_diagnostics) {
    # Stack diagnostics vertically
    p_mcmc_stacked <- (p_rhat / p_ess) +
        patchwork::plot_layout(heights = c(1, 1))

    grid_plots <- c(grid_plots, list(p_mcmc_stacked))
}


# Create the grid for standard plots
# Add the summary table to the end of the list
grid_plots <- c(grid_plots, list(p_table))

# Define layout design
# A-K are the 11 plots
# S is the stacked diagnostics (already in grid_plots as item 12)
# T is the summary table (item 13)
# 5 columns:
# Row 1: A B C D E
# Row 2: F G H I J
# Row 3: K S T T T (Table takes 3 slots)
layout_design <- "
ABCDE
FGHIJ
KSTTT
"

# Create the grid with custom design
p_grid <- patchwork::wrap_plots(grid_plots, design = layout_design)

# Right Panel: Settings (Diagnostics moved back to grid)
p_right <- p_settings

# Left Panel: Just the grid now (table is inside it)
p_left <- p_grid

# Combine: Left Panel | Right Panel
combined_plot <- (p_left | p_right) +
    patchwork::plot_layout(widths = c(6, 1)) +
    plot_annotation(title = plot_title, theme = plot_theme)

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
    width = 24, # Optimized for 24-inch screen
    height = 13.5, # 16:9 aspect ratio
    device = ragg::agg_png,
    bg = "white",
    dpi = 300 # High DPI for Retina display
)
