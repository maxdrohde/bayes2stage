# Generate correlated random intercept / slope longitudinal data

Generates simulated longitudinal data with correlated random intercepts
and slopes for testing two-stage design methods.

## Usage

``` r
generate_data(
  N = 2000,
  M = 5,
  alpha_main = 1,
  beta_x = 1,
  beta_z = 1,
  beta_t = 2,
  beta_t_x_interaction = 0.3,
  beta_t_z_interaction = 0,
  error_sd = 4,
  x_dist = c("normal", "poisson", "binomial", "negative_binomial", "beta_binomial"),
  x_size = NULL,
  x_disp_param = NULL,
  rand_intercept_sd = 3,
  rand_slope_sd = 1,
  rand_eff_corr = 0,
  gamma0 = 1,
  gamma1 = 1,
  gamma2 = 0,
  gamma_sd = 2
)
```

## Arguments

- N:

  Number of subjects (default: 2000)

- M:

  Number of time points per subject (default: 5)

- alpha_main:

  Intercept for the main model (default: 1)

- beta_x:

  Coefficient for x in the main model (default: 1)

- beta_z:

  Coefficient for z in the main model (default: 1)

- beta_t:

  Coefficient for time in the main model (default: 2)

- beta_t_x_interaction:

  Coefficient for time-x interaction (default: 0.3)

- beta_t_z_interaction:

  Coefficient for time-z interaction (default: 0)

- error_sd:

  Standard deviation of residual error (default: 4)

- x_dist:

  Distribution for x: "normal", "poisson", "binomial",
  "negative_binomial", or "beta_binomial" (default: "normal")

- x_size:

  Size parameter for binomial/beta_binomial distributions

- x_disp_param:

  Dispersion parameter for negative_binomial/beta_binomial

- rand_intercept_sd:

  Standard deviation of random intercepts (default: 3)

- rand_slope_sd:

  Standard deviation of random slopes (default: 1)

- rand_eff_corr:

  Correlation between random intercepts and slopes (default: 0)

- gamma0:

  Intercept for the imputation model (default: 1)

- gamma1:

  Linear coefficient for z in the imputation model (default: 1)

- gamma2:

  Quadratic coefficient for z in the imputation model (default: 0)

- gamma_sd:

  Standard deviation for the imputation model (normal distribution only)
  (default: 2)

## Value

A data frame with simulated longitudinal data
