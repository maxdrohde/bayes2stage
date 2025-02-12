#' Plot simulated data
#'
#' @param dataset Dataset to use
#' @param subset_size How many subjects to use?
#' @return A ggplot2 object
#' @export
data_plot <- function(dataset,
                      subset_size = 200){

  dataset <-
  dataset |>
  dplyr::filter(id %in% sample(unique(dataset$id), subset_size))

plt_t <-
  dataset |>
  ggplot2::ggplot() +
  ggplot2::aes(x = t, y = y, group = id) +
  ggplot2::geom_line(linewidth = 0.1) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Time", y = "Y")

plt_x <-
  dataset |>
  ggplot2::ggplot() +
  ggplot2::aes(x = x, y = y) +
  ggplot2::geom_point(size = 0.2) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "X", y = "Y")

plt_z <-
  dataset |>
  ggplot2::ggplot() +
  ggplot2::aes(x = z, y = y) +
  ggplot2::geom_point(size = 0.2) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Z", y = "Y")

plt_x <-
  dataset |>
  dplyr::filter(t == 1) |>
  ggplot2::ggplot() +
  ggplot2::aes(x = z, y = x) +
  ggplot2::geom_point(size = 0.2) +
  ggplot2::labs(x = "Z", y = "X")

hist_x <-
  dataset |>
  dplyr::filter(t == 1) |>
  ggplot2::ggplot() +
  ggplot2::aes(x) +
  ggplot2::geom_histogram(bins = 50, color = "black", fill = "lightgray") +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "X")

hist_z <-
  dataset |>
  dplyr::filter(t == 1) |>
  ggplot2::ggplot() +
  ggplot2::aes(z) +
  ggplot2::geom_histogram(bins = 50, color = "black", fill = "lightgray") +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Z")

design <- "AAAAAA
           BBCCDD
           EEEFFF"

combined <-
  patchwork::wrap_plots(
             A = plt_t,
             B = plt_x,
             C = plt_z,
             D = plt_x,
             E = hist_x,
             F = hist_z,
             design = design)

return(combined)
}
