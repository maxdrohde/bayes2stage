# Fit ACML model for outcome-dependent sampling

Fits an ascertainment-corrected maximum likelihood model for ODS
designs.

## Usage

``` r
fit_acml_ods(ods_df, cutoff_low, cutoff_high)
```

## Arguments

- ods_df:

  Data frame from ODS sampling

- cutoff_low:

  Lower cutoff value for sampling regions

- cutoff_high:

  Upper cutoff value for sampling regions

## Value

A data frame with parameter estimates and confidence intervals
