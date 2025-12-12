# Posterior Predictive Checks for Imputation Model

Functions to diagnose imputation model adequacy

## Usage

``` r
ppc_imputation(
  fit,
  data,
  imputation_covariates,
  distribution = c("normal", "bernoulli", "negative_binomial", "beta_binomial"),
  n_draws = 1000
)
```

## Arguments

- fit:

  CmdStanMCMC fit object from fit_stan_model()

- data:

  Data frame used for fitting

- imputation_covariates:

  Character vector of imputation model covariate names

- distribution:

  Imputation distribution used

- n_draws:

  Number of posterior draws for simulations

## Value

List with all diagnostic results and plots
