#' Given selected subjects, set X to missing for those not selected
#'
#' @param dataset Input data frame
#' @param selected_ids IDs for the selected subjects
#' @return data.frame
#' @export
set_missing <- function(dataset, selected_ids){

  # Check if required variables are present
  stopifnot("id must be a variable in the data" = ("id" %in% names(dataset)))
  stopifnot("x must be a variable in the data" = ("x" %in% names(dataset)))

  # Create indicator for being selected
  dataset$selected <- dataset$id %in% selected_ids

  # Set X values where not selected to NA
  dataset$x[!dataset$selected] <- NA

  return(dataset)
}
