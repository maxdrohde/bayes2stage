get_ods <- function(data,
                    sampling_type){

  # Convert to data.table
  data <- data.table::as.data.table(data)

  lms <-
    data[,
            .(lm_coef = list(stats::lm(y ~ t, data = .SD)$coefficients)),
            by = id]

  if (sampling_type == "intercept") {
    lms[, `:=`(target = lm_coef[[1]][["(Intercept)"]]), by = id]
    out <- lms[, .(id, target)]
  } else if (sampling_type == "slope") {
    lms[, `:=`(target = lm_coef[[1]][["t"]]), by = id]
    out <- lms[, .(id, target)]
  } else {
    stop("sampling_type must be 'intercept' or 'slope'")
  }

  out$sampling_type <- sampling_type

  return(tibble::as_tibble(out))
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
ods_design <- function(data,
                       sampling_type,
                       cutoff_high,
                       cutoff_low,
                       n_sampled,
                       prop_high,
                       prop_middle,
                       prop_low){

  stopifnot("Choose either 'intercept' or 'slope' as `sampling_type`" = sampling_type %in% c("intercept", "slope"))
  stopifnot("Strata proportions must sum to 1" = prop_high + prop_middle + prop_low == 1)
  stopifnot("x must be a variable in the data" = ("x" %in% names(data)))
  stopifnot("Number of subjects sampled must be a whole number" = is_whole_number(n_sampled))

  ods <- get_ods(data, sampling_type)

  # Categorize  into high, middle, and low strata based on quantiles
  ods$category <-
    cut(ods$target,
        breaks = c(-Inf,
                   stats::quantile(ods$target, cutoff_low),
                   stats::quantile(ods$target, cutoff_high),
                   Inf),
        labels = c("Low", "Middle", "High"))

  # Compute sampling sizes for each strata
  size_high <- n_sampled * prop_high
  size_middle <- n_sampled * prop_middle
  size_low <- n_sampled * prop_low

  # Check sizes
  sizes <- c(size_high, size_middle, size_low)
  stopifnot("Sample sizes must be positive whole numbers" =
              all(is_positive_integer(sizes)))

  # Convert sizes to integer
  size_high <- as.integer(size_high)
  size_middle <- as.integer(size_middle)
  size_low <- as.integer(size_low)

  # Sample IDs from the high, middle, and low strata
  high_ids <- sample(dplyr::filter(ods, category == "High")$id,
                     size = size_high,
                     replace = FALSE)
  middle_ids <- sample(dplyr::filter(ods, category == "Middle")$id,
                       size = size_middle,
                       replace = FALSE)
  low_ids <- sample(dplyr::filter(ods, category == "Low")$id,
                    size = size_low,
                    replace = FALSE)

  # Subjects IDs chosen for sampling
  selected_ids <- c(high_ids,
                    middle_ids,
                    low_ids)

  # If not chosen for stage 2, set x to missing
  stage2_df <- set_missing(data,
                           selected_ids)

  # Merge in the information on the estimated intercepts / slopes
  stage2_df <- dplyr::left_join(stage2_df,
                                ods,
                                by = "id")

  return(stage2_df)
}
