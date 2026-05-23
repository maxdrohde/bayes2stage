# Get default key parameters for diagnostics

Returns the Stan parameter names that are typically of interest for MCMC
diagnostics and model summaries. When `available` is provided, also
includes all indexed `beta[i]` and `gamma[i]` parameters found in the
available set (useful when models have multiple covariates).

## Usage

``` r
get_key_parameters(available = NULL)
```

## Arguments

- available:

  Optional character vector of available parameter names from a fitted
  model. When provided, all `beta[i]` and `gamma[i]` parameters are
  included.

## Value

Character vector of Stan parameter names
