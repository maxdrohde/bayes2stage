# Generate correlated random intercept / slope longitudinal data to test two-stage design methods. In the form of Schildcrout (2019) supplement

Generate correlated random intercept / slope longitudinal data to test
two-stage design methods. In the form of Schildcrout (2019) supplement

## Usage

``` r
generate_data_schild_2019_supp(
  N,
  Ms,
  alpha_main,
  beta_x,
  beta_z,
  beta_t,
  beta_t_x_interaction,
  beta_t_z_interaction,
  error_sd,
  x_prevalence,
  rand_intercept_sd,
  rand_slope_sd,
  rand_eff_corr,
  gamma0,
  gamma1,
  gamma2,
  gamma_sd
)
```

## Arguments

- N:

  Integer. Total number of subjects in the cohort. Default is 2000.

- Ms:

  Integer vector. Possible numbers of observations per subject, sampled
  uniformly. Default is c(4, 5, 6).

- alpha_main:

  Numeric. Fixed intercept (beta_0). Default is 75.

- beta_x:

  Numeric. Main effect of binary exposure X (beta_s). Default is -0.5.

- beta_z:

  Numeric. Effect of continuous confounder Z (beta_c). Default is -2.

- beta_t:

  Numeric. Fixed time effect (beta_t). Default is -1.

- beta_t_x_interaction:

  Numeric. Interaction between time and exposure (beta_st). Default is
  -0.5.

- beta_t_z_interaction:

  Numeric. Interaction between time and confounder. Default is 0.

- error_sd:

  Numeric. Standard deviation of measurement error (sigma_e). Default is
  sqrt(12.25) = 3.5.

- x_prevalence:

  Numeric. Prevalence of binary exposure X. Must be between 0 and 1.
  Default is 0.25.

- rand_intercept_sd:

  Numeric. Standard deviation of random intercept (sigma_0). Default is
  sqrt(81) = 9.

- rand_slope_sd:

  Numeric. Standard deviation of random slope (sigma_1). Default is
  sqrt(1.56) = 1.25.

- rand_eff_corr:

  Numeric. Correlation between random intercept and slope (rho). Must be
  between -1 and 1. Default is 0.

- gamma0:

  Numeric. Intercept for confounder Z generation. Default is 0.25.

- gamma1:

  Numeric. Linear coefficient for X in Z generation. Default is 0.5.

- gamma2:

  Numeric. Quadratic coefficient for X in Z generation. Default is 0.

- gamma_sd:

  Numeric. Standard deviation of Z residual error. Default is 1.

## Value

dataset
