#' Format the simulated data for Stan / NIMBLE
#'
#' @param dataset Dataset to use
#' @return A list suitable for input to MCMC software
#' @export
format_data_mcmc <- function(dataset,
                             spline = FALSE){

  # Check if required variables in dataset
  stopifnot("Some required variables not in the dataset" =
            all(c("y", "t", "id", "x", "z") %in% names(dataset)))

  id_df <-
    dataset |>
    dplyr::distinct(id, .keep_all = TRUE)

  X <- id_df[["x"]]
  Z <- id_df[["z"]]

  Xrep <- dataset[["x"]]
  Zrep <- dataset[["z"]]

  P <- ncol(X) + ncol(Z)
  S <- ncol(Z)

  data_list <-
    list(
      N = nrow(dataset),
      G = nrow(id_df),
      G_obs = sum(!is.na(id_df$x)),
      G_mis = sum(is.na(id_df$x)),
      X = X,
      Z = Z,
      Zrep = Zrep,
      P = P,
      S = S,
      t = dataset$t,
      y = dataset$y,
      index_obs = which(!is.na(id_df$x)),
      index_mis = which(is.na(id_df$x)),
      x_obs = id_df$x[!is.na(id_df$x)],
      id = dataset$id)

  return(data_list)
}
