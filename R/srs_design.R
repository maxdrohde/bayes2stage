#' Set x to missing based on an SRS design
#'
#' @param data Dataset to use
#' @param sampling_N How many subjects to sample?
#' @return A dataset where the x values are selected based on an SRS design
#' @export
srs_design <- function(data, sampling_N){

  # Check if arguments have correct type
  stopifnot("Number of subjects sampled must be integer" = is_positive_integer(sampling_N))

  # Check if required variables are present
  stopifnot("id must be a variable in the data" = ("id" %in% names(data)))
  stopifnot("x must be a variable in the data" = ("x" %in% names(data)))

  # IDs chosen for stage 2
  selected_ids <- sample(unique(data$id),
                         size = sampling_N,
                         replace = FALSE)

  # If not chosen for stage 2, set x to missing
  stage2_df <- set_missing(data,
                           selected_ids)

  return(stage2_df)
}
