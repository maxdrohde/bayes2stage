get_ods <- function(dataset, sampling_type){

  # Convert to data.table
  dataset <- data.table::as.data.table(dataset)

  lms <-
    dataset[,
            .(lm_coef = list(lm(y ~ t, data = .SD)$coefficients)),
            by = id]

  if (sampling_type == "intercept") {
    # Create an intercept column for each subject
    lms[, `:=`(intercept = lm_coef[[1]][["(Intercept)"]]), by = id]
    out <- lms[, .(id, intercept)]
  }

  if (sampling_type == "slope") {
    # Create a slope column for each subject
    lms[, `:=`(slope = lm_coef[[1]][["t"]]), by = id]
    out <- lms[, .(id, slope)]
  }

  return(tibble::as_tibble(out))
}

#' Set x to missing based on an ODS design
#'
#' @param dataset Dataset to use
#' @param sampling_type Which type of sampling? "intercept" or "slope"
#' @param cutoff_high Which quantile to use as the cutoff for the High category
#' @param cutoff_low Which quantile to use as the cutoff for the Low category
#' @param sampling_N How many subjects should be sampled?
#' @param prop_high What proportion to sample from the High category?
#' @param prop_middle What proportion to sample from the Middle category?
#' @param prop_low What proportion to sample from the Low category?
#' @return A dataset where the x values are selected based on an ODS design
#' @export
ods_design <- function(dataset,
                       sampling_type,
                       cutoff_high,
                       cutoff_low,
                       sampling_N,
                       prop_high,
                       prop_middle,
                       prop_low){

  stopifnot("Choose either 'intercept' or 'slope' as `sampling_type`" = sampling_type %in% c("intercept", "slope"))
  stopifnot("Strata proportions must sum to 1" = prop_high + prop_middle + prop_low == 1)
  stopifnot("x must be a variable in the data" = ("x" %in% names(dataset)))
  stopifnot("Number of subjects sampled must be a whole number" = is.wholenumber(sampling_N))

  ods <- get_ods(dataset, sampling_type)

  sampling_feature <- as.numeric(ods[,2])

  # Categorize  into high, middle, and low strata based on quantiles
  ods$category <-
    cut(sampling_feature,
        breaks = c(-Inf,
                   stats::quantile(sampling_feature, cutoff_low),
                   stats::quantile(sampling_feature, cutoff_high),
                   Inf),
        labels = c("Low", "Middle", "High"))

  size_high <- sampling_N * prop_high
  size_middle <- sampling_N * prop_middle
  size_low <- sampling_N * prop_low

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
  stage2_df <- dataset
  stage2_df[!(stage2_df$id %in% selected_ids),]$x <- NA

  # Mark which subjects were selected
  stage2_df$selected <- FALSE
  stage2_df[(stage2_df$id %in% selected_ids),]$selected <- TRUE

  # Merge in the information on the estimated intercepts / slopes
  stage2_df <- dplyr::left_join(stage2_df,
                                ods,
                                by = "id")

  return(stage2_df)
}
