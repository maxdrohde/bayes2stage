get_ods_old <- function(dataset, type){

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
