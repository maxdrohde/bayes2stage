# Parallel coordinates plot for diagnosing divergences

Parallel coordinates plot for diagnosing divergences

## Usage

``` r
plot_parcoord(fit, parameters = NULL)
```

## Arguments

- fit:

  CmdStanMCMC fit object

- parameters:

  Character vector of parameter names. If NULL, uses parameters related
  to random effects.

## Value

A ggplot object
