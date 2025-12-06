# Extract MCMC summary statistics

A wrapper around MCMCvis::MCMCsummary to extract summary statistics from
MCMC samples.

## Usage

``` r
mcmc_summary(mcmc_output, dataset_id)
```

## Arguments

- mcmc_output:

  MCMC output object (e.g., from fit_model)

- dataset_id:

  A character string identifying the dataset (added as a column)

## Value

A data frame with summary statistics for each parameter
