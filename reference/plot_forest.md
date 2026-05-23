# Forest plot with table-style annotations

Creates a forest plot with point estimates, confidence intervals, and
right-aligned numeric columns for SE and CI text, using a `theme_void()`
layout with a manual axis.

## Usage

``` r
plot_forest(
  data,
  ref_line = 0,
  xlab = NULL,
  panel_label = NULL,
  ci_format = "%5.2f (%5.2f, %4.2f)",
  se_format = "%4.2f",
  extra_col = NULL,
  font_family = "sans",
  mono_family = "mono"
)
```

## Arguments

- data:

  Data frame with required columns: `label`, `est`, `lo`, `hi`, `se`.
  Row types are inferred automatically:

  - **Data row**: `est` is not `NA` – gets point + CI + text columns.

  - **Header row**: `est` is `NA` but `label` is non-empty – bold,
    flush-left, no CI.

  - **Spacer row**: `label` is `""` or `NA` – empty vertical gap.

  Data rows appearing before the first header are rendered flush-left;
  data rows after a header are indented.

- ref_line:

  X-intercept for vertical reference line, or `NULL` to omit.

- xlab:

  X-axis label, or `NULL` to omit.

- panel_label:

  Panel annotation (e.g., `"A)"`) placed top-left, or `NULL` to omit.

- ci_format:

  [`sprintf()`](https://rdrr.io/r/base/sprintf.html) format string
  applied to `(est, lo, hi)` for the CI text column.

- se_format:

  [`sprintf()`](https://rdrr.io/r/base/sprintf.html) format string
  applied to `se` for the SE text column.

- extra_col:

  Optional named list specifying an additional numeric column placed
  between the SE and CI columns, or `NULL` to omit. Required elements:
  `col` (name of the column in `data`), `format`
  ([`sprintf()`](https://rdrr.io/r/base/sprintf.html) format string),
  and `header` (column header text).

- font_family:

  Font family for labels, headers, and axis text.

- mono_family:

  Monospace font family for numeric columns.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object.
