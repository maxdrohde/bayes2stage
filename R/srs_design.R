#' Set x to missing based on an SRS design
#'
#' @param dataset Dataset to use
#' @param sampling_N What number of subjects to sample?
#' @return A dataset where the x values are selected based on an SRS design
#' @export
srs_design <- function(dataset, sampling_N){

  stopifnot("id must be a variable in the data" = ("id" %in% names(dataset)))
  stopifnot("x must be a variable in the data" = ("x" %in% names(dataset)))
  stopifnot("Number of subjects sampled must be a whole number" = is.wholenumber(sampling_N))


  # IDs chosen for stage 2
  selected_ids <- sample(unique(dataset$id),
                         size = sampling_N,
                         replace = FALSE)

  # If not chosen for stage 2, set x to missing
  stage2_df <- dataset
  stage2_df[!(stage2_df$id %in% selected_ids),]$x <- NA

  # Mark which subjects were selected
  stage2_df$selected <- TRUE
  stage2_df[(stage2_df$id %in% selected_ids),]$selected <- FALSE

  return(stage2_df)
}
