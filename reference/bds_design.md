# Set x to missing based on an BDS design

Set x to missing based on an BDS design

## Usage

``` r
bds_design(
  data,
  fixed_effects_formula,
  sampling_type,
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

- fixed_effects_formula:

  Formula for the fixed-effects when fitting the model to estimate BLUPs

- sampling_type:

  Which type of sampling? "intercept" or "slope"

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

A dataset where the x values are selected based on an BDS design
