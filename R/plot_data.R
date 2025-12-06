#' Plot simulated data
#'
#' @param data Dataset to use
#' @param subset_size How many subjects to use?
#' @return A ggplot2 object
#' @export
plot_data <- function(data,
                      subset_size = 200){

  data <-
  data |>
  dplyr::filter(id %in% sample(unique(data$id), subset_size))

plt_t <-
  data |>
  ggplot2::ggplot() +
  ggplot2::aes(x = t, y = y, group = id) +
  ggplot2::geom_line(linewidth = 0.1) +
  ggplot2::geom_hline(yintercept = 0, linetype = 2) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Time", y = "Y")

plt_xy <-
  data |>
  ggplot2::ggplot() +
  ggplot2::aes(x = x, y = y) +
  ggplot2::geom_point(size = 0.2) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "X", y = "Y") +
  ggplot2::coord_equal()

plt_zy <-
  data |>
  ggplot2::ggplot() +
  ggplot2::aes(x = z, y = y) +
  ggplot2::geom_point(size = 0.2) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Z", y = "Y") +
  ggplot2::coord_equal()

plt_zx <-
  data |>
  dplyr::filter(t == 1) |>
  ggplot2::ggplot() +
  ggplot2::aes(x = z, y = x) +
  ggplot2::geom_point(size = 0.2) +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "Z", y = "X") +
  ggplot2::coord_equal()

hist_x <-
  data |>
  dplyr::filter(t == 1) |>
  ggplot2::ggplot() +
  ggplot2::aes(x) +
  ggplot2::geom_histogram(bins = 50, color = "black", fill = "lightgray") +
  ggplot2::theme_bw() +
  ggplot2::labs(x = "X")

hist_z <-
  data |>
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
             B = plt_xy,
             C = plt_zy,
             D = plt_zx,
             E = hist_x,
             F = hist_z,
             design = design)

return(combined)
}
