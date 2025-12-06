# Plot MCMC trace plots

A wrapper around MCMCvis::MCMCtrace to create trace plots of MCMC
samples.

## Usage

``` r
mcmc_trace(mcmc_output, print_to_pdf = FALSE)
```

## Arguments

- mcmc_output:

  MCMC output object (e.g., from fit_model)

- print_to_pdf:

  Logical; if TRUE, saves trace plots to a PDF file

## Value

Trace plots of MCMC samples
