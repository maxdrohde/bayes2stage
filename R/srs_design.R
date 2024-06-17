#' Set x_e to missing based on an SRS design
#'
#' @param dataset Dataset to use
#' @param proportion What proportion of subjects to sample?
#' @return A dataset where the x_e values are selected based on an SRS design
#' @export
srs_design <- function(dataset, proportion){

  ##############################################################################
  # Check if proper names in dataset
  col_names <- names(dataset)
  stopifnot("y" %in% col_names)
  stopifnot("t" %in% col_names)
  stopifnot("id" %in% col_names)
  stopifnot("x_e" %in% col_names)
  ##############################################################################

  # IDs chosen for stage 2
  selected_ids <- sample(unique(dataset$id),
                         size = proportion * length(unique(dataset$id)),
                         replace = FALSE)

  # If not chosen for stage 2, set x_e to missing
  stage2_df <- dataset
  stage2_df[!(stage2_df$id %in% selected_ids),]$x_e <- NA

  # Mark selected = TRUE if chosen for stage 2 collection
  stage2_df$selected <- "No"
  stage2_df[(stage2_df$id %in% selected_ids),]$selected <- "Yes"

  return(stage2_df)
}
