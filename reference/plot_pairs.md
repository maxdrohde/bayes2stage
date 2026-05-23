# Pairs plots for diagnosing posterior geometry

Pairs plots for diagnosing posterior geometry

## Usage

``` r
plot_pairs(fit, parameters = NULL)
```

## Arguments

- fit:

  CmdStanMCMC fit object

- parameters:

  Character vector of parameter names to plot. If NULL, uses a subset of
  key parameters.

## Value

A ggplot object
