# PPC for count imputation model (negative binomial)

PPC for count imputation model (negative binomial)

## Usage

``` r
ppc_imputation_count(fit, data, imputation_covariates, n_draws = 500)
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

List with rootogram and dispersion check
