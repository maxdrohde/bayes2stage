# Changelog

## bayes2stage 0.0.0.9000

- Stan-based model fitting via
  [`fit_stan_model()`](https://maxdrohde.github.io/bayes2stage/reference/fit_stan_model.md)
  using the ‘instantiate’ package.
- Support for four covariate distributions: normal, bernoulli, negative
  binomial, and beta binomial.
- Three sampling design functions:
  [`srs_design()`](https://maxdrohde.github.io/bayes2stage/reference/srs_design.md),
  [`ods_design()`](https://maxdrohde.github.io/bayes2stage/reference/ods_design.md),
  [`bds_design()`](https://maxdrohde.github.io/bayes2stage/reference/bds_design.md).
- Ascertainment-corrected maximum likelihood estimation via
  [`fit_acml_ods()`](https://maxdrohde.github.io/bayes2stage/reference/fit_acml_ods.md).
- Data simulation with
  [`generate_data()`](https://maxdrohde.github.io/bayes2stage/reference/generate_data.md)
  supporting multiple covariate distributions.
- Added package vignette with complete workflow example.
- Legacy NIMBLE-based fitting available via
  [`fit_model()`](https://maxdrohde.github.io/bayes2stage/reference/fit_model.md).
