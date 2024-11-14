# Given a dataset and an lmer formula, return a data.frame of the BLUPs for
# each subject ID
get_blups <- function(dataset, type){


  mod <- lme4::lmer("y ~ x_z + t + (1 + t | id)",
                    data = dataset)

  if (type == "intercept") {
    out <-
      lme4::ranef(mod) |>
      as.data.frame() |>
      dplyr::filter(term == "(Intercept)") |>
      dplyr::select(id = grp, intercept = condval)
  }

  if (type == "slope") {
    out <-
      lme4::ranef(mod) |>
      as.data.frame() |>
      dplyr::filter(term == "t") |>
      dplyr::select(id = grp, slope = condval)
  }

  out$id <-
    out$id |>
    as.character() |>
    as.integer()

  return(as.data.frame(out))
}

#' Set x_e to missing based on an BDS design
#'
#' @param dataset Dataset to use
#' @param type Which type of sampling? "intercept" or "slope"
#' @param cutoff_high Which quantile to use as the cutoff for the High category
#' @param cutoff_low Which quantile to use as the cutoff for the Low category
#' @param prop_high What proportion to sample from the High category?
#' @param prop_middle What proportion to sample from the Middle category?
#' @param prop_low What proportion to sample from the Low category?
#' @return A dataset where the x_e values are selected based on an BDS design
#' @export
bds_design <- function(dataset,
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

  bds <- get_blups(dataset, type)
  G <- nrow(bds)

  target <- as.numeric(bds[,2])

  # Categorize  into high, middle, and low strata based on quantiles
  bds$category <-
    cut(target,
        breaks = c(-Inf,
                   stats::quantile(target, cutoff_low),
                   stats::quantile(target, cutoff_high),
                   Inf),
        labels = c("Low", "Middle", "High"))

  # Sample IDs from the high, middle, and low strata
  high_ids <- sample(dplyr::filter(bds, category == "High")$id,
                     size = G * prop_high,
                     replace = FALSE)
  middle_ids <- sample(dplyr::filter(bds, category == "Middle")$id,
                       size = G * prop_middle,
                       replace = FALSE)
  low_ids <- sample(dplyr::filter(bds, category == "Low")$id,
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
                                bds,
                                by = "id")

  return(stage2_df)
}
