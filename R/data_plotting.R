#' Plot simulated data
#'
#' @param dataset Dataset to use
#' @param subset_size How many subjects to use?
#' @return A ggplot2 object
#' @export
data_plot <- function(dataset, subset_size){

  df <-
  df |>
  dplyr::filter(id %in% sample(unique(df$id), subset_size))

plt_t <-
  df |>
  ggplot2::ggplot() +
  ggplot2::aes(x = t, y = y, group = id) +
  ggplot2::geom_line(linewidth = 0.1) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::theme_bw()

plt_xe <-
  df |>
  ggplot2::ggplot() +
  ggplot2::aes(x = x_e, y = y) +
  ggplot2::geom_point(size = 0.2) +
  ggplot2::theme_bw()

plt_xz <-
  df |>
  ggplot2::ggplot() +
  ggplot2::aes(x = x_z, y = y) +
  ggplot2::geom_point(size = 0.2) +
  ggplot2::theme_bw()

plt_x <-
  df |>
  dplyr::filter(t == 1) |>
  ggplot2::ggplot() +
  ggplot2::aes(x = x_z, y = x_e) +
  ggplot2::geom_point(size = 0.2) +
  ggplot2::theme_bw()

design <- "AAA
           BCD"

combined <-
  patchwork::wrap_plots(
             A = plt_t,
             B = plt_xe,
             C = plt_xz,
             D = plt_x,
             design = design)

return(combined)
}
