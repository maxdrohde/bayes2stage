# Format the simulated data for Stan

Format the simulated data for Stan

## Usage

``` r
format_data_mcmc(
  data,
  main_model_formula = NULL,
  imputation_model_formula = NULL,
  imputation_distribution = c("normal", "bernoulli", "beta_binomial",
    "negative_binomial")
)
```

## Arguments

- data:

  Dataset to use

- main_model_formula:

  One-sided formula or string for covariates in the main model (e.g.,
  `~ age + splines::ns(bmi, 3)`). Intercept is automatically removed.

- imputation_model_formula:

  One-sided formula or string for covariates in the imputation model
  (e.g., `~ age + factor(site)`). Intercept is automatically removed.

- imputation_distribution:

  Distribution for the imputation model: "normal" for continuous x,
  "bernoulli" for binary x, "beta_binomial" for bounded count data, or
  "negative_binomial" for unbounded count data

## Value

A list suitable for input to MCMC software
