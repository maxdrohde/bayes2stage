# Given a dataset and an `lmer` formula, return a data.frame
# of the BLUPs for each subject
get_blups <- function(data,
                      fixed_effects_formula,
                      sampling_type){

  check_cols(data, c("t", "id"))

  if (!("t" %in% all.vars(fixed_effects_formula))) {
    cli::cli_abort("{.var t} must be in {.arg fixed_effects_formula}.")
  }

  # Add the random-effects to the user supplied fixed-effects formula
  lmer_formula <-
    glue::glue("{deparse(fixed_effects_formula)} + (1 + t | id)") |>
    stats::as.formula()

  # Fit the lme4 mixed-effect model to compute BLUPs
  mod <- lme4::lmer(formula = lmer_formula,
                    data = data,
                    REML = FALSE)

  # Extact the BLUPs
  if (sampling_type == "intercept") {
    out <-
      lme4::ranef(mod) |>
      as.data.frame() |>
      dplyr::filter(term == "(Intercept)") |>
      dplyr::select(id = grp,
                    target = condval)
  } else if (sampling_type == "slope") {
    out <-
      lme4::ranef(mod) |>
      as.data.frame() |>
      dplyr::filter(term == "t") |>
      dplyr::select(id = grp,
                    target = condval)
  } else {
    cli::cli_abort("{.arg sampling_type} must be {.val intercept} or {.val slope}.")
  }

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
bds_design <- function(data,
                       fixed_effects_formula,
                       sampling_type,
                       cutoff_high,
                       cutoff_low,
                       n_sampled,
                       prop_high,
                       prop_middle,
                       prop_low){

  if (!(sampling_type %in% c("intercept", "slope"))) {
    cli::cli_abort("{.arg sampling_type} must be {.val intercept} or {.val slope}.")
  }
  if (prop_high + prop_middle + prop_low != 1) {
    cli::cli_abort("Strata proportions must sum to 1.")
  }
  check_cols(data, "x")
  if (!is_positive_integer(n_sampled)) {
    cli::cli_abort("{.arg n_sampled} must be a positive integer.")
  }

  blups <-
    get_blups(data,
              fixed_effects_formula,
              sampling_type)

  # Categorize  into high, middle, and low strata based on quantiles
  blups$category <-
    cut(blups$target,
        breaks = c(-Inf,
                   stats::quantile(blups$target, cutoff_low),
                   stats::quantile(blups$target, cutoff_high),
                   Inf),
        labels = c("Low", "Middle", "High"))

  # Compute sampling sizes for each strata
  size_high <- n_sampled * prop_high
  size_middle <- n_sampled * prop_middle
  size_low <- n_sampled * prop_low

  # Check sizes
  sizes <- c(size_high, size_middle, size_low)
  if (!all(is_positive_integer(sizes))) {
    cli::cli_abort("Sample sizes must be positive whole numbers.")
  }

  # Convert sizes to integer
  size_high <- as.integer(size_high)
  size_middle <- as.integer(size_middle)
  size_low <- as.integer(size_low)

  # Sample IDs from the high, middle, and low strata
  high_ids <- sample(dplyr::filter(blups, category == "High")$id,
                     size = size_high,
                     replace = FALSE)

  middle_ids <- sample(dplyr::filter(blups, category == "Middle")$id,
                       size = size_middle,
                       replace = FALSE)

  low_ids <- sample(dplyr::filter(blups, category == "Low")$id,
                    size = size_low,
                    replace = FALSE)

  # Subjects IDs chosen for sampling
  selected_ids <- c(high_ids,
                    middle_ids,
                    low_ids)

  # If not chosen for stage 2, set x to missing
  stage2_df <- set_missing(data,
                           selected_ids)

  # Merge in the information on the estimated BLUPs
  stage2_df <- dplyr::left_join(stage2_df,
                                blups,
                                by = "id")

  return(stage2_df)
}
