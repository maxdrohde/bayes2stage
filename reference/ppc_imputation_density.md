# Posterior predictive density check for imputation model

Posterior predictive density check for imputation model

## Usage

``` r
ppc_imputation_density(fit, data, imputation_covariates, n_draws = 100)
```

## Arguments

- fit:

  CmdStanMCMC fit object from fit_stan_model()

- data:

  Data frame used for fitting

- imputation_covariates:

  Character vector of imputation model covariates

- n_draws:

  Number of posterior draws to visualize

## Value

ggplot object
