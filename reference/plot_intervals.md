# Posterior interval plots

Posterior interval plots

## Usage

``` r
plot_intervals(fit, parameters = NULL, prob = 0.5, prob_outer = 0.95)
```

## Arguments

- fit:

  CmdStanMCMC fit object

- parameters:

  Character vector of parameter names to plot. If NULL, uses default key
  parameters.

- prob:

  Inner probability interval (default: 0.5)

- prob_outer:

  Outer probability interval (default: 0.95)

## Value

A ggplot object
