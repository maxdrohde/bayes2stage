#' Perform stratified sampling based on a target variable,
#' such as subject-specific intercepts or slopes
#'
#' @param id_target_df Data frame with `id` and `target` columns
#' @param n_sampled Total number of subjects to sample
#' @param cutoff_low Quantile for Low category cutoff
#' @param cutoff_high Quantile for High category cutoff
#' @param prop_low Proportion to sample from Low category
#' @param prop_middle Proportion to sample from Middle category
#' @param prop_high Proportion to sample from High category
#' @return A list with two elements:
#'   - `selected_ids`: Vector of selected subject IDs
#'   - `id_target_df`: The input data frame with `category` column added
#' @noRd
perform_stratified_sampling <- function(
    id_target_df,
    n_sampled,
    cutoff_low,
    cutoff_high,
    prop_low,
    prop_middle,
    prop_high
) {
    # Validate inputs
    if (abs(prop_high + prop_middle + prop_low - 1) > .Machine$double.eps^0.5) {
        cli::cli_abort("Strata proportions must sum to 1.")
    }
    if (!is_positive_integer(n_sampled)) {
        cli::cli_abort("{.arg n_sampled} must be a positive integer.")
    }

    # Categorize into high, middle, and low strata based on quantiles
    breaks <- c(
        -Inf,
        stats::quantile(id_target_df$target, cutoff_low),
        stats::quantile(id_target_df$target, cutoff_high),
        Inf
    )

    id_target_df$category <- cut(
        id_target_df$target,
        breaks = breaks,
        labels = c("Low", "Middle", "High")
    )

    # Compute sampling sizes for each strata
    sizes <- c(
        High = as.integer(n_sampled * prop_high),
        Middle = as.integer(n_sampled * prop_middle),
        Low = as.integer(n_sampled * prop_low)
    )

    if (!all(is_positive_integer(sizes))) {
        cli::cli_abort("Sample sizes must be positive whole numbers.")
    }

    # Validate that each stratum has enough candidates
    strata_counts <- table(id_target_df$category)
    for (stratum in names(sizes)) {
        available <- as.integer(strata_counts[stratum])
        if (is.na(available)) available <- 0L
        requested <- sizes[stratum]
        if (available < requested) {
            cli::cli_abort(
                "Stratum {.val {stratum}} has {available} subjects but {requested} were requested."
            )
        }
    }

    # Helper to sample IDs from a specific category
    sample_ids <- function(cat, size) {
        candidates <- dplyr::filter(id_target_df, category == cat)$id
        sample(candidates, size = size, replace = FALSE)
    }

    selected_ids <- c(
        sample_ids("High", sizes["High"]),
        sample_ids("Middle", sizes["Middle"]),
        sample_ids("Low", sizes["Low"])
    )

    result <- list(
        selected_ids = selected_ids,
        id_target_df = id_target_df
    )

    return(result)
}
