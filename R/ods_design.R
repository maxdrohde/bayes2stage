get_ods <- function(dataset, type){

  # Nest by ID
  # then fit lm() to each row
  lms <-
    dataset |>
    dplyr::nest_by(id) |>
    dplyr::summarise(lm_output = list(stats::lm("y ~ t", data = data)),
                     .groups = "drop")

  if (type == "intercept") {
    # Extract the intercept from each fitted lm() object
    lms$intercept <-
      purrr::map_dbl(lms$lm_output,
                     ~.x$coefficients[["(Intercept)"]])

    out <-
      lms |>
      dplyr::select(id, intercept)
  }

  if (type == "slope") {
    # Extract the slope from each fitted lm() object
    lms$slope <-
      purrr::map_dbl(lms$lm_output,
                     ~.x$coefficients[["t"]])

    out <-
      lms |>
      dplyr::select(id, slope)
  }

  return(as.data.frame(out))
}


#' Set x_e to missing based on an ODS design
#'
#' @param dataset Dataset to use
#' @param type Which type of sampling? "intercept" or "slope"
#' @param cutoff_high Which quantile to use as the cutoff for the High category
#' @param cutoff_low Which quantile to use as the cutoff for the Low category
#' @param prop_high What proportion to sample from the High category?
#' @param prop_middle What proportion to sample from the Middle category?
#' @param prop_low What proportion to sample from the Low category?
#' @return A dataset where the x_e values are selected based on an ODS design
#' @export
ods_design <- function(dataset,
                       type,
                       cutoff_high = 0.9,
                       cutoff_low = 0.1,
                       prop_high = 0.1,
                       prop_middle = 0.05,
                       prop_low = 0.1){

  ##############################################################################
  # Check if proper names in dataset
  col_names <- names(dataset)
  stopifnot("y" %in% col_names)
  stopifnot("t" %in% col_names)
  stopifnot("id" %in% col_names)
  stopifnot("x_e" %in% col_names)

  stopifnot("Choose either 'intercept' or 'slope' as `type`" = type %in% c("intercept", "slope"))
  ##############################################################################

  ods <- get_ods(dataset, type)
  G <- nrow(ods)

  target <- as.numeric(ods[,2])

  # Categorize  into high, middle, and low strata based on quantiles
  ods$category <-
    cut(target,
        breaks = c(-Inf,
                   stats::quantile(target, 0.1),
                   stats::quantile(target, 0.9),
                   Inf),
        labels = c("Low", "Middle", "High"))

  # Sample IDs from the high, middle, and low strata
  high_ids <- sample(dplyr::filter(ods, category == "High")$id,
                     size = G * prop_high,
                     replace = FALSE)
  middle_ids <- sample(dplyr::filter(ods, category == "Middle")$id,
                       size = G * prop_middle,
                       replace = FALSE)
  low_ids <- sample(dplyr::filter(ods, category == "Low")$id,
                    size = G * prop_low,
                    replace = FALSE)

  # IDs chosen for stage 2
  selected_ids <- c(high_ids, middle_ids, low_ids)

  # If not chosen for stage 2, set x_e to missing
  stage2_df <- dataset
  stage2_df[!(stage2_df$id %in% selected_ids),]$x_e <- NA

  # Mark selected = TRUE if chosen for stage 2 collection
  stage2_df$selected <- "No"
  stage2_df[(stage2_df$id %in% selected_ids),]$selected <- "Yes"

  stage2_df <- dplyr::left_join(stage2_df,
                                ods,
                                by = "id")

  return(stage2_df)
}
