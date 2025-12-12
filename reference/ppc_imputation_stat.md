# Posterior predictive check for a summary statistic

Posterior predictive check for a summary statistic

## Usage

``` r
ppc_imputation_stat(
  fit,
  data,
  imputation_covariates,
  stat = mean,
  stat_name = "Mean",
  n_draws = 1000
)
```

## Arguments

- fit:

  CmdStanMCMC fit object

- data:

  Data frame used for fitting

- imputation_covariates:

  Character vector of imputation model covariates

- stat:

  Function to compute statistic

- stat_name:

  Name for display

- n_draws:

  Number of posterior draws

## Value

List with observed value, predicted distribution, p-value, and plot
