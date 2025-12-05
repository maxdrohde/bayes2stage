#' Format the simulated data for Stan / NIMBLE
#'
#' @param dataset Dataset to use
#' @return A list suitable for input to MCMC software
#' @export
format_data_mcmc <- function(dataset,
                             main_vars = NULL,   # columns in the main model
                             imputation_vars = NULL)   # columns in the imputation model
  {

  dataset <-
    dataset |>
    dplyr::arrange(id, t)

  # Re-index id to 1..G
  dataset <-
    dataset |>
    dplyr::mutate(id_idx = as.integer(factor(id)))

  id_df <-
    dataset |>
    dplyr::distinct(id_idx,
                    .keep_all = TRUE) |>
    dplyr::arrange(id_idx)

  X <- as.matrix(dataset[, main_vars, drop = FALSE])
  Z <- as.matrix(id_df[, imputation_vars, drop = FALSE])

  P <- ncol(X)
  S <- ncol(Z)

  data_list <-
    list(
      N       = nrow(dataset),
      G       = nrow(id_df),
      G_obs   = sum(!is.na(id_df$x)),
      G_mis   = sum(is.na(id_df$x)),
      P       = P,
      S       = S,
      t       = dataset$t,
      X       = X,
      Z       = Z,
      y       = dataset$y,
      index_obs = which(!is.na(id_df$x)),
      index_mis = which(is.na(id_df$x)),
      x_obs     = id_df$x[!is.na(id_df$x)],
      id        = dataset$id_idx
    )

  return(data_list)
}
