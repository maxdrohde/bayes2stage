# PPC for Bernoulli imputation model

PPC for Bernoulli imputation model

## Usage

``` r
ppc_imputation_bernoulli(fit, data, imputation_covariates, n_draws = 1000)
```

## Arguments

- fit:

  CmdStanMCMC fit object

- data:

  Data frame

- imputation_covariates:

  Character vector of imputation model covariates

- n_draws:

  Number of posterior draws

## Value

List with calibration plot, ROC metrics
