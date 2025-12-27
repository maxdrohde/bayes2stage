get_ods <- function(data, sampling_type) {
  validate_sampling_type(sampling_type)

  coef_name <- c(intercept = "(Intercept)", slope = "t")[[sampling_type]]

  out <- data |>
    dplyr::summarize(
      target = stats::lm(y ~ t)$coefficients[[coef_name]],
      .by = id
    )
  out$sampling_type <- sampling_type

  return(out)
}

#' Set x to missing based on an ODS design
#'
#' @param data Dataset to use
#' @param sampling_type Which type of sampling? "intercept" or "slope"
#' @param cutoff_high Which quantile to use as the cutoff for the High category
#' @param cutoff_low Which quantile to use as the cutoff for the Low category
#' @param n_sampled How many subjects should be sampled?
#' @param prop_high What proportion to sample from the High category?
#' @param prop_middle What proportion to sample from the Middle category?
#' @param prop_low What proportion to sample from the Low category?
#' @return A dataset where the x values are selected based on an ODS design
#' @export
ods_design <- function(
  data,
  sampling_type,
  cutoff_high,
  cutoff_low,
  n_sampled,
  prop_high,
  prop_middle,
  prop_low
) {
  check_cols(data, "x")

  ods <- get_ods(data, sampling_type)

  sampling_result <- perform_stratified_sampling(
    id_target_df = ods,
    n_sampled = n_sampled,
    cutoff_low = cutoff_low,
    cutoff_high = cutoff_high,
    prop_low = prop_low,
    prop_middle = prop_middle,
    prop_high = prop_high
  )

  # If not chosen for stage 2, set x to missing
  stage2_df <- set_missing(data, sampling_result$selected_ids)

  # Merge in the information on the estimated intercepts / slopes (includes category)
  stage2_df <- dplyr::left_join(stage2_df, sampling_result$id_target_df, by = "id")

  return(stage2_df)
}
