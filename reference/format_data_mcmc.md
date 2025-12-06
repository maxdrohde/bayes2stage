# Format the simulated data for Stan / NIMBLE

Format the simulated data for Stan / NIMBLE

## Usage

``` r
format_data_mcmc(
  data,
  main_model_covariates = NULL,
  imputation_model_covariates = NULL,
  imputation_distribution = c("normal", "bernoulli", "beta_binomial",
    "negative_binomial")
)
```

## Arguments

- data:

  Dataset to use

- main_model_covariates:

  Character vector of column names for covariates in the main model

- imputation_model_covariates:

  Character vector of column names for covariates in the imputation
  model

- imputation_distribution:

  Distribution for the imputation model: "normal" for continuous x,
  "bernoulli" for binary x, "beta_binomial" for bounded count data, or
  "negative_binomial" for unbounded count data

## Value

A list suitable for input to MCMC software
