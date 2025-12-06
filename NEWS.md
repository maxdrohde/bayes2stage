# bayes2stage 0.0.0.9000

* Stan-based model fitting via `fit_stan_model()` using the 'instantiate' package.
* Support for four covariate distributions: normal, bernoulli, negative binomial, and beta binomial.
* Three sampling design functions: `srs_design()`, `ods_design()`, `bds_design()`.
* Ascertainment-corrected maximum likelihood estimation via `fit_acml_ods()`.
* Data simulation with `generate_data()` supporting multiple covariate distributions.
* Added package vignette with complete workflow example.
* Legacy NIMBLE-based fitting available via `fit_model()`.
