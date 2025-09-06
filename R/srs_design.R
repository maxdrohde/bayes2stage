#' Set x to missing based on an SRS design
#'
#' @param dataset Dataset to use
#' @param N How many subjects to sample?
#' @return A dataset where the x values are selected based on an SRS design
#' @export
srs_design <- function(dataset, N){

  # Check if arguments have correct type
  stopifnot("Number of subjects sampled must be integer" = is_positive_integer(N))

  # Check if required variables are present
  stopifnot("id must be a variable in the data" = ("id" %in% names(dataset)))
  stopifnot("x must be a variable in the data" = ("x" %in% names(dataset)))

  # IDs chosen for stage 2
  selected_ids <- sample(unique(dataset$id),
                         size = N,
                         replace = FALSE)

  # If not chosen for stage 2, set x to missing
  stage2_df <- set_missing(dataset, selected_ids)

  return(stage2_df)
}
