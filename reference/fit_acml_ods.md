# Fit ACML model for outcome-dependent sampling

Fits an ascertainment-corrected maximum likelihood model for ODS
designs.

## Usage

``` r
fit_acml_ods(
  ods_df,
  cutoff_low,
  cutoff_high,
  main_model_formula = ~1,
  random_effects_formula = ~1 + t
)
```

## Arguments

- ods_df:

  Data frame from ODS sampling. Must contain columns: `id`, `target`,
  `selected`, `category`, `x`, `y`, `t`, `sampling_type`.

- cutoff_low:

  Lower quantile cutoff (between 0 and 1) for sampling regions.

- cutoff_high:

  Upper quantile cutoff (between 0 and 1) for sampling regions.

- main_model_formula:

  One-sided formula or string for additional covariates in the main
  model (e.g., `~ z`). Terms `t`, `x`, and `x:t` are always included.
  Default: `~ 1` (no additional covariates).

- random_effects_formula:

  Formula for random effects. Default: `~ 1 + t` (random intercept and
  slope on time).

## Value

A data frame with class `"acml_fit"` containing parameter estimates and
confidence intervals. Columns: `variable`, `mean`, `sd`, `2.5%`,
`97.5%`.
