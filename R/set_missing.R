#' Given selected subjects, set X to missing for those not selected
#'
#' @param data Input data frame
#' @param selected_ids IDs for the selected subjects
#' @return data.frame
#' @export
set_missing <- function(data, selected_ids){

  check_cols(data, c("id", "x"))

  # Create indicator for being selected
  data$selected <- data$id %in% selected_ids

  # Set X values where not selected to NA
  data$x[!data$selected] <- NA

  return(data)
}
