# Package index

## All functions

- [`ACi1q()`](https://maxdrohde.github.io/bayes2stage/reference/ACi1q.md)
  : Ascertainment correction piece for univariate sampling

- [`ACi2q()`](https://maxdrohde.github.io/bayes2stage/reference/ACi2q.md)
  : Ascertainment correction piece for bivariate sampling

- [`CreateSubjectData()`](https://maxdrohde.github.io/bayes2stage/reference/CreateSubjectData.md)
  : Create a list of subject-specific data

- [`LogLikeC.Score2()`](https://maxdrohde.github.io/bayes2stage/reference/LogLikeC.Score2.md)
  : Calculate the gradient of the conditional likelihood for the
  univariate and bivariate sampling cases across all subjects
  (CheeseCalc=FALSE) or the cheese part of the sandwich estimator if
  CheeseCalc=TRUE.

- [`LogLikeC2()`](https://maxdrohde.github.io/bayes2stage/reference/LogLikeC2.md)
  : Calculate the conditional likelihood for the univariate and
  bivariate sampling cases across all subjects (Keep.liC=FALSE) or the
  subject specific contributions to the conditional likelihood along
  with the log-transformed ascertainment correction for multiple
  imputation (Keep.liC=TRUE).

- [`LogLikeCAndScore2()`](https://maxdrohde.github.io/bayes2stage/reference/LogLikeCAndScore2.md)
  : Calculate the ascertainment corrected log likelihood and score

- [`LogLikeiC2()`](https://maxdrohde.github.io/bayes2stage/reference/LogLikeiC2.md)
  : Calculate the ss contributions to the conditional likelihood for the
  univariate and bivariate sampling cases.

- [`acml.lmem2()`](https://maxdrohde.github.io/bayes2stage/reference/acml.lmem2.md)
  : Fitting function: ACML or WL for a linear mixed effects model
  (random intercept and slope)

- [`bds_design()`](https://maxdrohde.github.io/bayes2stage/reference/bds_design.md)
  : Set x to missing based on an BDS design

- [`cantor_seed()`](https://maxdrohde.github.io/bayes2stage/reference/cantor_seed.md)
  : Generate a unique seed from two indices using Cantor pairing

- [`check_cols()`](https://maxdrohde.github.io/bayes2stage/reference/check_cols.md)
  : Check for required columns in a data frame

- [`fit_acml_ods()`](https://maxdrohde.github.io/bayes2stage/reference/fit_acml_ods.md)
  : Fit ACML model for outcome-dependent sampling

- [`fit_stan_model()`](https://maxdrohde.github.io/bayes2stage/reference/fit_stan_model.md)
  : Fit a Bayesian two-stage model using Stan

- [`format_data_mcmc()`](https://maxdrohde.github.io/bayes2stage/reference/format_data_mcmc.md)
  : Format the simulated data for Stan

- [`generate_data()`](https://maxdrohde.github.io/bayes2stage/reference/generate_data.md)
  : Generate correlated random intercept / slope longitudinal data

- [`generate_data_schild_2019_supp()`](https://maxdrohde.github.io/bayes2stage/reference/generate_data_schild_2019_supp.md)
  : Generate correlated random intercept / slope longitudinal data to
  test two-stage design methods. In the form of Schildcrout (2019)
  supplement

- [`li.lme()`](https://maxdrohde.github.io/bayes2stage/reference/li.lme.md)
  : Calculate a subject-specific contribution to a log-likelihood for
  longitudinal normal data

- [`li.lme.score2()`](https://maxdrohde.github.io/bayes2stage/reference/li.lme.score2.md)
  : Subject specific contribution to the lme model score (also returns
  marginal Vi=Cov(Y\|X))

- [`logACi1q()`](https://maxdrohde.github.io/bayes2stage/reference/logACi1q.md)
  : Log of the Ascertainment correction for univariate sampling

- [`logACi1q.score2()`](https://maxdrohde.github.io/bayes2stage/reference/logACi1q.score2.md)
  : Gradient of the log of the ascertainment correction piece for
  sampling based on univariate Q_i

- [`logACi2q()`](https://maxdrohde.github.io/bayes2stage/reference/logACi2q.md)
  : Log of the Ascertainment correction piece for bivariate sampling

- [`logACi2q.score2()`](https://maxdrohde.github.io/bayes2stage/reference/logACi2q.score2.md)
  : Gradient of the log of the ascertainment correction piece for
  sampling based on bivariate Q_i.

- [`mcmc_forest()`](https://maxdrohde.github.io/bayes2stage/reference/mcmc_forest.md)
  : Create a forest plot of MCMC output

- [`mcmc_summary()`](https://maxdrohde.github.io/bayes2stage/reference/mcmc_summary.md)
  : Extract MCMC summary statistics

- [`mcmc_trace()`](https://maxdrohde.github.io/bayes2stage/reference/mcmc_trace.md)
  : Plot MCMC trace plots

- [`ods_design()`](https://maxdrohde.github.io/bayes2stage/reference/ods_design.md)
  : Set x to missing based on an ODS design

- [`plot_data()`](https://maxdrohde.github.io/bayes2stage/reference/plot_data.md)
  : Plot simulated data

- [`ppc_imputation()`](https://maxdrohde.github.io/bayes2stage/reference/ppc_imputation.md)
  : Posterior Predictive Checks for Imputation Model

- [`ppc_imputation_bernoulli()`](https://maxdrohde.github.io/bayes2stage/reference/ppc_imputation_bernoulli.md)
  : PPC for Bernoulli imputation model

- [`ppc_imputation_beta_binomial()`](https://maxdrohde.github.io/bayes2stage/reference/ppc_imputation_beta_binomial.md)
  : PPC for beta-binomial imputation model

- [`ppc_imputation_count()`](https://maxdrohde.github.io/bayes2stage/reference/ppc_imputation_count.md)
  : PPC for count imputation model (negative binomial)

- [`ppc_imputation_density()`](https://maxdrohde.github.io/bayes2stage/reference/ppc_imputation_density.md)
  : Posterior predictive density check for imputation model

- [`ppc_imputation_intervals()`](https://maxdrohde.github.io/bayes2stage/reference/ppc_imputation_intervals.md)
  : Check prediction interval coverage

- [`ppc_imputation_residuals()`](https://maxdrohde.github.io/bayes2stage/reference/ppc_imputation_residuals.md)
  : Residual diagnostics for normal imputation model

- [`ppc_imputation_stat()`](https://maxdrohde.github.io/bayes2stage/reference/ppc_imputation_stat.md)
  : Posterior predictive check for a summary statistic

- [`print(`*`<ppc_imputation>`*`)`](https://maxdrohde.github.io/bayes2stage/reference/print.ppc_imputation.md)
  : Print method for ppc_imputation results

- [`set_missing()`](https://maxdrohde.github.io/bayes2stage/reference/set_missing.md)
  : Given selected subjects, set X to missing for those not selected

- [`srs_design()`](https://maxdrohde.github.io/bayes2stage/reference/srs_design.md)
  : Set x to missing based on an SRS design

- [`vi.calc()`](https://maxdrohde.github.io/bayes2stage/reference/vi.calc.md)
  :

  Calculate V_i = Z_i D t(Z_i) + sig_e^2 `I_(n_i)`
