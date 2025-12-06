
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayes2stage

<!-- badges: start -->

<!-- badges: end -->

`bayes2stage` provides tools for designing and analyzing Bayesian
two-stage studies with longitudinal data. The package supports
outcome-dependent sampling (ODS) and BLUP-dependent sampling (BDS)
designs, where expensive covariates are measured on a subset of subjects
selected based on their outcomes or random effects. Bayesian analysis is
performed using Stan for efficient mixed effects modeling with
imputation.

## Key Features

- **Simulate longitudinal data** with random intercepts/slopes and
  various covariate distributions
- **Implement two-stage sampling designs**: Simple random sampling
  (SRS), outcome-dependent sampling (ODS), and BLUP-dependent sampling
  (BDS)
- **Bayesian model fitting** using Stan with automatic imputation for
  missing covariates
- **Support for multiple distributions**: Normal and Bernoulli
  covariates in the imputation model
- **Flexible model specification** for both main outcome and imputation
  models

## Installation

You can install the development version of bayes2stage from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("maxdrohde/bayes2stage")
```

## Example Workflow

Here’s a complete example demonstrating the two-stage design workflow:

``` r
library(bayes2stage)

# Step 1: Generate longitudinal data
# Create a dataset with N=2000 subjects, M=5 time points per subject
# x is a continuous covariate that will be expensive to measure
set.seed(123)
full_data <- generate_data(
  N = 2000,
  M = 5,
  alpha_main = 1,
  beta_x = 1,
  beta_z = 1,
  beta_t = 2,
  beta_t_x_interaction = 0.3,
  x_dist = "normal",
  rand_intercept_sd = 3,
  rand_slope_sd = 1,
  rand_eff_corr = 0.5
)

# Step 2: Apply an outcome-dependent sampling design
# Select 200 subjects for stage 2 based on their estimated random slopes
# Sample more heavily from the tails of the slope distribution
stage2_data <- ods_design(
  dataset = full_data,
  sampling_type = "slope",
  cutoff_high = 0.75,
  cutoff_low = 0.25,
  n_sampled = 200,
  prop_high = 0.4,
  prop_middle = 0.2,
  prop_low = 0.4
)

# Step 3: Fit Bayesian mixed effects model with imputation
# The model will impute missing x values for subjects not selected in stage 2
fit <- fit_stan_model(
  data = stage2_data,
  main_model_covariates = c("x", "z", "t"),
  imputation_model_covariates = c("z"),
  imputation_distribution = "normal",
  n_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000
)

# Step 4: Examine results
fit$summary(variables = c("beta_x", "beta_z", "beta_t"))
```

## Two-Stage Sampling Designs

The package implements three sampling designs for selecting which
subjects to measure the expensive covariate on:

- **SRS** (`srs_design()`): Simple random sampling - select subjects at
  random
- **ODS** (`ods_design()`): Outcome-dependent sampling - select based on
  subject-specific intercepts or slopes estimated from simple linear
  regressions
- **BDS** (`bds_design()`): BLUP-dependent sampling - select based on
  best linear unbiased predictors (BLUPs) from a mixed effects model

All designs support stratified sampling, allowing you to oversample from
the tails of the distribution while maintaining some subjects from the
middle.

## Model Fitting

The package supports two approaches for Bayesian model fitting:

- **Stan-based** (`fit_stan_model()`): Modern, efficient HMC sampling
  using Stan via the `instantiate` package (recommended)
- **NIMBLE-based** (`fit_model()`): Alternative implementation using
  NIMBLE (legacy support)

The Stan models handle missing covariate data through a joint model that
simultaneously: 1. Models the outcome as a function of covariates and
random effects 2. Imputes missing covariate values based on observed
covariates

## Citation

If you use this package in your research, please cite:

    Rohde, M. (2024). bayes2stage: Bayesian Analysis of Two-Stage Designs.
    R package version 0.0.0.9000. https://github.com/maxdrohde/bayes2stage

## Getting Help

- Report bugs or request features:
  <https://github.com/maxdrohde/bayes2stage/issues>
- Questions about usage: Open a discussion on GitHub
