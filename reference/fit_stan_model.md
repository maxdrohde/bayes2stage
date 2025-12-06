# Fit a Bayesian two-stage model using Stan

Fits a mixed effects model with imputation using Stan via the
instantiate package.

## Usage

``` r
fit_stan_model(
  data,
  main_model_covariates,
  imputation_model_covariates,
  imputation_distribution = c("normal", "bernoulli", "beta_binomial",
    "negative_binomial"),
  n_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.8,
  seed = 777L,
  parallel_chains = 1L
)
```

## Arguments

- data:

  A data frame containing the outcome and covariates

- main_model_covariates:

  Character vector of covariate names for the main model

- imputation_model_covariates:

  Character vector of covariate names for the imputation model

- imputation_distribution:

  Distribution for the imputation model: "normal" for continuous x,
  "bernoulli" for binary x, "beta_binomial" for bounded count data, or
  "negative_binomial" for unbounded count data (default: "normal")

- n_chains:

  Number of MCMC chains (default: 4)

- iter_warmup:

  Number of warmup iterations per chain (default: 1000)

- iter_sampling:

  Number of sampling iterations per chain (default: 1000)

- adapt_delta:

  Target acceptance rate for HMC (default: 0.8)

- seed:

  Random seed for reproducibility (default: 777L)

- parallel_chains:

  Number of chains to run in parallel (default: 1L)

## Value

A CmdStanMCMC fit object
