# Set x to missing based on an ellipse design

Stratifies subjects by Mahalanobis distance in the bivariate (intercept,
slope) space. Subjects far from the center (high Mahalanobis distance)
are oversampled.

## Usage

``` r
ellipse_design(
  data,
  estimation_method = c("ods", "bds"),
  fixed_effects_formula = NULL,
  cutoff_low,
  cutoff_high,
  n_sampled,
  prop_high,
  prop_middle,
  prop_low
)
```

## Arguments

- data:

  Dataset to use

- estimation_method:

  How to estimate per-subject intercepts and slopes: `"ods"` uses
  per-subject OLS, `"bds"` uses BLUPs from a mixed model.

- fixed_effects_formula:

  Formula for the fixed-effects when fitting the model to estimate
  BLUPs. Required when `estimation_method = "bds"`.

- cutoff_low:

  Which quantile to use as the cutoff for the Low category

- cutoff_high:

  Which quantile to use as the cutoff for the High category

- n_sampled:

  How many subjects should be sampled?

- prop_high:

  What proportion to sample from the High category?

- prop_middle:

  What proportion to sample from the Middle category?

- prop_low:

  What proportion to sample from the Low category?

## Value

A dataset where the x values are selected based on an ellipse design.
Includes columns `target` (Mahalanobis distance), `category`,
`intercept`, and `slope`.
