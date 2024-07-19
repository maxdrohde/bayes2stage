#' Format the simulated data for Stan
#'
#' @param dataset Dataset to use
#' @return A list suitable for input to Stan
#' @export
create_stan_data <- function(dataset,
                             spline = FALSE){

  ##############################################################################
  # Check if proper names in dataset
  col_names <- names(dataset)
  stopifnot("y" %in% col_names)
  stopifnot("t" %in% col_names)
  stopifnot("id" %in% col_names)
  stopifnot("x_e" %in% col_names)
  ##############################################################################

  id_df <-
    dataset |>
    dplyr::filter(t == 1)

  X <- dataset[c("x_z")] |> as.matrix()
  t <- dataset$t
  P <- ncol(X)

  if (spline) {
    Z <- model.matrix(~-1 + splines::ns(x_z, 4), data = id_df)
  } else {
    Z <- id_df[c("x_z")] |> as.matrix()
  }

  S <- ncol(Z)

  data_list <-
    list(
      N = nrow(dataset),
      G = nrow(id_df),
      G_obs = sum(!is.na(id_df$x_e)),
      G_mis = sum(is.na(id_df$x_e)),
      P = P,
      X = X,
      t = t,
      S = S,
      Z = Z,
      y = dataset$y,
      index_obs = which(!is.na(id_df$x_e)),
      index_mis = which(is.na(id_df$x_e)),
      x_e_obs = id_df$x_e[!is.na(id_df$x_e)],
      id = dataset$id)

  return(data_list)
}
