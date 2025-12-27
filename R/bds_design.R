# Given a dataset and an `lmer` formula, return a data.frame
# of the BLUPs for each subject
get_blups <- function(data, fixed_effects_formula, sampling_type) {
  check_cols(data, c("t", "id"))
  validate_sampling_type(sampling_type)

  if (!("t" %in% all.vars(fixed_effects_formula))) {
    cli::cli_abort("{.var t} must be in {.arg fixed_effects_formula}.")
  }

  lmer_formula <-
    glue::glue("{deparse(fixed_effects_formula)} + (1 + t | id)") |>
    stats::as.formula()

  mod <- lme4::lmer(formula = lmer_formula, data = data, REML = FALSE)

  term_to_filter <- c(intercept = "(Intercept)", slope = "t")[[sampling_type]]

  out <-
    lme4::ranef(mod) |>
    as.data.frame() |>
    dplyr::filter(term == term_to_filter) |>
    dplyr::select(id = grp, target = condval)

  out$sampling_type <- sampling_type
  out$id <- as.integer(as.character(out$id))

  return(out)
}

#' Set x to missing based on an BDS design
#'
#' @param data Dataset to use
#' @param fixed_effects_formula Formula for the fixed-effects when fitting the model to estimate BLUPs
#' @param sampling_type Which type of sampling? "intercept" or "slope"
#' @param cutoff_high Which quantile to use as the cutoff for the High category
#' @param cutoff_low Which quantile to use as the cutoff for the Low category
#' @param n_sampled How many subjects should be sampled?
#' @param prop_high What proportion to sample from the High category?
#' @param prop_middle What proportion to sample from the Middle category?
#' @param prop_low What proportion to sample from the Low category?
#' @return A dataset where the x values are selected based on an BDS design
#' @export
bds_design <- function(
  data,
  fixed_effects_formula,
  sampling_type,
  cutoff_high,
  cutoff_low,
  n_sampled,
  prop_high,
  prop_middle,
  prop_low
) {
  check_cols(data, "x")

  blups <-
    get_blups(data, fixed_effects_formula, sampling_type)

  sampling_result <- perform_stratified_sampling(
    id_target_df = blups,
    n_sampled = n_sampled,
    cutoff_low = cutoff_low,
    cutoff_high = cutoff_high,
    prop_low = prop_low,
    prop_middle = prop_middle,
    prop_high = prop_high
  )

  # If not chosen for stage 2, set x to missing
  stage2_df <- set_missing(data, sampling_result$selected_ids)

  # Merge in the information on the estimated BLUPs (includes category)
  stage2_df <- dplyr::left_join(stage2_df, sampling_result$id_target_df, by = "id")

  return(stage2_df)
}
