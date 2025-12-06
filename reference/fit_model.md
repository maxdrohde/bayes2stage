# Fit a Bayesian two-stage model using NIMBLE

Fits a mixed effects model with imputation using NIMBLE for MCMC
sampling.

## Usage

``` r
fit_model(
  data,
  main_model_covariates,
  imputation_model_covariates,
  imputation_model_distribution,
  correlated_random_effects = TRUE,
  n_chains = 4,
  niter = 10000,
  nburnin = 2000,
  x_size = NULL,
  print_summary = FALSE,
  print_code = FALSE
)
```

## Arguments

- data:

  A data frame containing the outcome and covariates. Must include
  columns: y (outcome), t (time), x (exposure), and id (subject
  identifier).

- main_model_covariates:

  Character vector of covariate names for the main model

- imputation_model_covariates:

  Character vector of covariate names for the imputation model

- imputation_model_distribution:

  Distribution for the imputation model. One of: "normal", "binomial",
  "beta_binomial", "poisson", "negative_binomial"

- correlated_random_effects:

  Logical; if TRUE, random intercepts and slopes are correlated
  (default: TRUE)

- n_chains:

  Number of MCMC chains (default: 4)

- niter:

  Number of MCMC iterations per chain (default: 10000)

- nburnin:

  Number of burn-in iterations to discard (default: 2000)

- x_size:

  Integer vector of trial sizes; required for binomial or beta_binomial
  distributions (default: NULL)

- print_summary:

  Logical; if TRUE, prints MCMC summary (default: FALSE)

- print_code:

  Logical; if TRUE, prints the NIMBLE model code (default: FALSE)

## Value

A list containing MCMC samples and optionally WAIC
