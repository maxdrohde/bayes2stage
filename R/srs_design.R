#' Set x to missing based on an SRS design
#'
#' @param data Dataset to use
#' @param n_sampled How many subjects to sample?
#' @return A dataset where the x values are selected based on an SRS design
#' @export
srs_design <- function(data, n_sampled){

  if (!is_positive_integer(n_sampled)) {
    cli::cli_abort("{.arg n_sampled} must be a positive integer.")
  }
  check_cols(data, c("id", "x"))

  n_subjects <- length(unique(data$id))
  if (n_sampled > n_subjects) {
    cli::cli_abort("{.arg n_sampled} ({n_sampled}) exceeds number of subjects ({n_subjects}).")
  }

  # IDs chosen for stage 2
  all_ids <- unique(data$id)
  selected_ids <- all_ids[sample.int(length(all_ids), size = n_sampled)]

  # If not chosen for stage 2, set x to missing
  stage2_df <- set_missing(data,
                           selected_ids)

  return(stage2_df)
}
