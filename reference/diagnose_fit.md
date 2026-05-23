# Comprehensive MCMC diagnostics

Runs a full suite of MCMC diagnostics and returns plots and numeric
summaries.

## Usage

``` r
diagnose_fit(fit, parameters = NULL, print_summary = TRUE)
```

## Arguments

- fit:

  CmdStanMCMC fit object from fit_stan_model()

- parameters:

  Character vector of parameter names to diagnose. If NULL, uses default
  key parameters.

- print_summary:

  Logical; print diagnostic summary to console (default: TRUE)

## Value

A list with plots and diagnostic summaries
