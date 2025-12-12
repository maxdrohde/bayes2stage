# PPC for beta-binomial imputation model

PPC for beta-binomial imputation model

## Usage

``` r
ppc_imputation_beta_binomial(fit, data, imputation_covariates, n_draws = 500)
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

List with diagnostics
