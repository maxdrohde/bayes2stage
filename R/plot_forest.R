#' Forest plot with table-style annotations
#'
#' Creates a forest plot with point estimates, confidence intervals, and
#' right-aligned numeric columns for SE and CI text, using a `theme_void()`
#' layout with a manual axis.
#'
#' @param data Data frame with required columns: `label`, `est`, `lo`, `hi`,
#'   `se`. Row types are inferred automatically:
#'   - **Data row**: `est` is not `NA` -- gets point + CI + text columns.
#'   - **Header row**: `est` is `NA` but `label` is non-empty -- bold,
#'     flush-left, no CI.
#'   - **Spacer row**: `label` is `""` or `NA` -- empty vertical gap.
#'
#'   Data rows appearing before the first header are rendered flush-left;
#'   data rows after a header are indented.
#' @param ref_line X-intercept for vertical reference line, or `NULL` to omit.
#' @param xlab X-axis label, or `NULL` to omit.
#' @param panel_label Panel annotation (e.g., `"A)"`) placed top-left, or
#'   `NULL` to omit.
#' @param ci_format [sprintf()] format string applied to `(est, lo, hi)` for
#'   the CI text column.
#' @param se_format [sprintf()] format string applied to `se` for the SE text
#'   column.
#' @param extra_col Optional named list specifying an additional numeric column
#'   placed between the SE and CI columns, or `NULL` to omit. Required elements:
#'   `col` (name of the column in `data`), `format` ([sprintf()] format string),
#'   and `header` (column header text).
#' @param font_family Font family for labels, headers, and axis text.
#' @param mono_family Monospace font family for numeric columns.
#'
#' @return A [ggplot2::ggplot] object.
#'
#' @export
plot_forest <- function(
    data,
    ref_line = 0,
    xlab = NULL,
    panel_label = NULL,
    ci_format = "%5.2f (%5.2f, %4.2f)",
    se_format = "%4.2f",
    extra_col = NULL,
    font_family = "sans",
    mono_family = "mono"
) {
    required_cols <- c("label", "est", "lo", "hi", "se")
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0L) {
        cli::cli_abort("Missing required column{?s}: {.field {missing_cols}}")
    }

    n <- nrow(data)

    # Classify row types
    is_header <- is.na(data$est) & !is.na(data$label) & data$label != ""
    is_spacer <- is.na(data$label) | data$label == ""
    is_data   <- !is.na(data$est)

    # y positions (top row = n, bottom row = 1)
    data$y <- n:1L

    # Formatted text columns
    data$se_str <- ifelse(is_data, sprintf(se_format, data$se), "")
    data$ci_str <- ifelse(
        is_data,
        sprintf(ci_format, data$est, data$lo, data$hi),
        ""
    )

    # --- Axis layout (proportional to axis range) ---
    data_range <- range(c(data$lo, data$hi), na.rm = TRUE)
    breaks     <- pretty(data_range)
    axis_min   <- min(breaks)
    axis_max   <- max(breaks)
    axis_range <- axis_max - axis_min

    has_extra <- !is.null(extra_col)
    if (has_extra) {
        extra_values <- data[[extra_col$col]]
        data$extra_str <- ifelse(
            is_data & !is.na(extra_values),
            sprintf(extra_col$format, extra_values),
            ""
        )
    }

    header_x   <- axis_min - 0.5625 * axis_range
    indent_x   <- axis_min - 0.4375 * axis_range
    se_x       <- axis_max + 0.3125 * axis_range
    extra_x    <- axis_max + 0.5625 * axis_range
    ci_x       <- axis_max + (if (has_extra) 1.125 else 0.875) * axis_range
    xlim_left  <- axis_min - 0.75   * axis_range
    xlim_right <- axis_max + (if (has_extra) 1.625 else 1.375) * axis_range

    # Label x-position: rows before the first header are flush-left
    first_header <- which(is_header)[1L]
    if (is.na(first_header)) {
        data$x_lab <- header_x
    } else {
        data$x_lab <- ifelse(
            seq_len(n) < first_header | is_header,
            header_x,
            indent_x
        )
    }

    # --- Subsets ---
    dat_ci      <- data[is_data, ]
    dat_headers <- data[is_header, ]
    dat_plain   <- data[is_data, ]

    # --- y reference points ---
    y_min        <- min(data$y)
    y_max        <- max(data$y)
    axis_y       <- y_min - 1.2
    tick_bottom  <- y_min - 1.7
    tick_label_y <- y_min - 2.3
    title_y      <- y_min - 3.5
    col_header_y <- y_max + 1.5
    underline_y  <- y_max + 1.0

    # --- Build plot ---
    p <- ggplot2::ggplot(data, ggplot2::aes(y = y))

    # Reference line
    if (!is.null(ref_line)) {
        p <- p + ggplot2::annotate(
            "segment",
            x = ref_line, xend = ref_line,
            y = axis_y, yend = y_max + 0.5,
            linewidth = 0.3
        )
    }

    p <- p +
        # CI error bars
        ggplot2::geom_errorbar(
            data = dat_ci,
            ggplot2::aes(xmin = lo, xmax = hi),
            width = 0.3, linewidth = 0.4
        ) +
        # Point estimates
        ggplot2::geom_point(
            data = dat_ci,
            ggplot2::aes(x = est),
            size = 1.5
        )

    # Labels: data rows (plain)
    if (nrow(dat_plain) > 0L) {
        p <- p + ggplot2::geom_text(
            data = dat_plain,
            ggplot2::aes(x = x_lab, label = label),
            hjust = 0, size = 3, family = font_family
        )
    }

    # Labels: header rows (bold)
    if (nrow(dat_headers) > 0L) {
        p <- p + ggplot2::geom_text(
            data = dat_headers,
            ggplot2::aes(x = x_lab, label = label),
            hjust = 0, size = 3, fontface = "bold", family = font_family
        )
    }

    p <- p +
        # SE column
        ggplot2::geom_text(
            data = dat_ci,
            ggplot2::aes(label = se_str),
            x = se_x, hjust = 0.5, size = 2.8, family = mono_family
        ) +
        # CI column
        ggplot2::geom_text(
            data = dat_ci,
            ggplot2::aes(label = ci_str),
            x = ci_x, hjust = 0.5, size = 2.8, family = mono_family
        )

    # Extra column (between SE and CI)
    if (has_extra) {
        p <- p +
            ggplot2::geom_text(
                data = dat_ci,
                ggplot2::aes(label = extra_str),
                x = extra_x, hjust = 0.5, size = 2.8, family = mono_family
            )
    }

    # --- Column headers ---
    p <- p +
        ggplot2::annotate(
            "text", x = se_x, y = col_header_y, label = "SE",
            hjust = 0.5, fontface = "bold", size = 3.5, family = font_family
        ) +
        ggplot2::annotate(
            "text", x = ci_x, y = col_header_y, label = "Estimate (95% CI)",
            hjust = 0.5, fontface = "bold", size = 3.5, family = font_family
        ) +
        # Header underlines
        ggplot2::annotate(
            "segment",
            x    = se_x - 0.125 * axis_range,
            xend = se_x + 0.125 * axis_range,
            y = underline_y, yend = underline_y, linewidth = 0.4
        ) +
        ggplot2::annotate(
            "segment",
            x    = ci_x - 0.375 * axis_range,
            xend = ci_x + 0.375 * axis_range,
            y = underline_y, yend = underline_y, linewidth = 0.4
        )

    if (has_extra) {
        p <- p +
            ggplot2::annotate(
                "text", x = extra_x, y = col_header_y, label = extra_col$header,
                hjust = 0.5, fontface = "bold", size = 3.5, family = font_family
            ) +
            ggplot2::annotate(
                "segment",
                x    = extra_x - 0.125 * axis_range,
                xend = extra_x + 0.125 * axis_range,
                y = underline_y, yend = underline_y, linewidth = 0.4
            )
    }

    # --- Panel label ---
    if (!is.null(panel_label)) {
        p <- p + ggplot2::annotate(
            "text",
            x = xlim_left + 0.0625 * axis_range, y = y_max + 2.5,
            label = panel_label, hjust = 0, size = 4, family = font_family
        )
    }

    # --- Custom x-axis ---
    p <- p +
        # Axis line
        ggplot2::annotate(
            "segment",
            x = axis_min, xend = axis_max,
            y = axis_y, yend = axis_y, linewidth = 0.5
        ) +
        # Tick marks
        ggplot2::annotate(
            "segment",
            x = breaks, xend = breaks,
            y = axis_y, yend = tick_bottom, linewidth = 0.5
        ) +
        # Tick labels
        ggplot2::annotate(
            "text", x = breaks, y = tick_label_y,
            label = as.character(breaks), size = 3, family = font_family
        )

    # Axis title
    if (!is.null(xlab)) {
        p <- p + ggplot2::annotate(
            "text",
            x = (axis_min + axis_max) / 2, y = title_y,
            label = xlab, size = 3.5, family = font_family
        )
    }

    # --- Scales and theme ---
    ylim_bottom <- if (!is.null(xlab)) title_y - 1 else tick_label_y - 1
    ylim_top    <- if (!is.null(panel_label)) y_max + 3.0 else col_header_y + 1.0

    p <- p +
        ggplot2::scale_x_continuous(
            limits = c(xlim_left, xlim_right), expand = c(0, 0)
        ) +
        ggplot2::scale_y_continuous(
            limits = c(ylim_bottom, ylim_top), expand = c(0, 0)
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.background = ggplot2::element_rect(fill = "white", color = NA)
        )

    p
}
