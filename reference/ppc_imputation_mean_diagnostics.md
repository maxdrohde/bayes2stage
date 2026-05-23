# Mean model diagnostics for normal imputation

Comprehensive diagnostics to detect misspecification of the conditional
mean E\[x\|z\] in the normal imputation model. Produces four diagnostic
plots that help identify nonlinearity, systematic patterns, and model
fit issues.

## Usage

``` r
ppc_imputation_mean_diagnostics(fit, data, imputation_covariates, n_bins = 10)
```

## Arguments

- fit:

  CmdStanMCMC fit object from fit_stan_model()

- data:

  Data frame used for fitting

- imputation_covariates:

  Character vector of imputation model covariate names

- n_bins:

  Number of bins for binned residual plot (default: 10)

## Value

A list containing:

- plots:

  List of ggplot objects: resid_vs_z, partial_residuals, obs_vs_pred,
  binned_residuals

- r_squared:

  R-squared of observed vs predicted

- fitted:

  Vector of fitted values E\[x\|z\]

- residuals:

  Vector of residuals (observed - fitted)
