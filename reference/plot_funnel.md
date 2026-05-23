# Quick funnel diagnostic for hierarchical parameters

Creates a scatter plot of log(sigma) vs random effects to diagnose
funnel geometry.

## Usage

``` r
plot_funnel(
  fit,
  sigma_param = "sigma_re[1]",
  effect_param = "re\\[1,",
  n_effects = 5
)
```

## Arguments

- fit:

  CmdStanMCMC fit object

- sigma_param:

  Name of the sigma parameter (default: `"sigma_re[1]"`)

- effect_param:

  Pattern for random effect parameters (default: `"re\\[1,"`)

- n_effects:

  Number of random effects to plot (default: 5)

## Value

A ggplot object
