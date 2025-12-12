# Check prediction interval coverage

Check prediction interval coverage

## Usage

``` r
ppc_imputation_intervals(
  fit,
  data,
  imputation_covariates,
  prob = 0.95,
  n_draws = 1000
)
```

## Arguments

- fit:

  CmdStanMCMC fit object

- data:

  Data frame

- imputation_covariates:

  Character vector of imputation model covariates

- prob:

  Coverage probability

- n_draws:

  Number of posterior draws

## Value

List with coverage proportion and plot
