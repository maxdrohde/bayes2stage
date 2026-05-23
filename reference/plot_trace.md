# Trace plots for MCMC diagnostics

Trace plots for MCMC diagnostics

## Usage

``` r
plot_trace(fit, parameters = NULL)
```

## Arguments

- fit:

  CmdStanMCMC fit object from fit_stan_model()

- parameters:

  Character vector of parameter names to plot. If NULL, uses default key
  parameters.

## Value

A ggplot object
