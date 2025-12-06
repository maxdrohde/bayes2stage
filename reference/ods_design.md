# Set x to missing based on an ODS design

Set x to missing based on an ODS design

## Usage

``` r
ods_design(
  data,
  sampling_type,
  cutoff_high,
  cutoff_low,
  n_sampled,
  prop_high,
  prop_middle,
  prop_low
)
```

## Arguments

- data:

  Dataset to use

- sampling_type:

  Which type of sampling? "intercept" or "slope"

- cutoff_high:

  Which quantile to use as the cutoff for the High category

- cutoff_low:

  Which quantile to use as the cutoff for the Low category

- n_sampled:

  How many subjects should be sampled?

- prop_high:

  What proportion to sample from the High category?

- prop_middle:

  What proportion to sample from the Middle category?

- prop_low:

  What proportion to sample from the Low category?

## Value

A dataset where the x values are selected based on an ODS design
