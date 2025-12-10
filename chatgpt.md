# Paper Outline: Bayesian Two-Stage Designs for Longitudinal Studies with Costly Covariates

## 1. Aim and Positioning

- Goal: demonstrate efficient inference when an expensive, subject-level
  covariate is measured only on a subsample, leveraging longitudinal
  outcomes and Bayesian joint modeling.
- Audience: methodological statisticians and applied researchers
  designing biomarker/genomics/clinical monitoring studies.
- Deliverables: evidence that outcome/BLUP-dependent designs + joint
  Bayesian imputation improve efficiency vs SRS and frequentist ACML.

## 2. Introduction & Motivation

- Problem: high-cost covariates; abundant repeated outcomes; need to
  maximize information under cost constraints.
- Two-stage designs: collect y,z on everyone; measure x on selected
  subjects; oversample trajectory tails (ODS/BDS) to gain power for x
  effects.
- Gaps: limited Bayesian tools that (a) handle multiple x
  distributions, (b) fully marginalize discrete x, (c) integrate design
  selection, (d) provide open-source implementations.
- Contributions of `bayes2stage`:
  - Unified Stan models for normal, Bernoulli, negative binomial,
    beta-binomial x with random intercept/slope outcome model and x
    imputation.
  - Optimizations (pos/len indices, marginalized discrete x, vectorized
    likelihoods) for computational efficiency.
  - Design functions:
    [`srs_design()`](https://maxdrohde.github.io/bayes2stage/reference/srs_design.md),
    [`ods_design()`](https://maxdrohde.github.io/bayes2stage/reference/ods_design.md)
    (OLS intercept/slope strata),
    [`bds_design()`](https://maxdrohde.github.io/bayes2stage/reference/bds_design.md)
    (BLUP strata).
  - Frequentist comparator: ascertainment-corrected ML
    ([`fit_acml_ods()`](https://maxdrohde.github.io/bayes2stage/reference/fit_acml_ods.md)).
  - Open-source R package + vignette for full workflow.

## 3. Methods

- Data structure & outcome model: subjects i=1..N, times j=1..M; y model
  with fixed effects (x, z, t, x\*t), random intercept/slope with
  correlation; t scaled 0–1; z observed for all; x partially observed.
- Stage-2 selection strategies:
  - SRS baseline.
  - ODS: strata on OLS intercept/slope; cutoff_low/high and stratum
    proportions.
  - BDS: strata on BLUP intercept/slope from mixed model; sensitivity to
    M and RE variance.
- Bayesian joint model:
  - Outcome: non-centered random effects, correlated intercept/slope;
    priors exponential on scales, weakly-informative normals on fixed
    effects.
  - Imputation families: normal (continuous x), Bernoulli (logit),
    negative binomial (log-mean with dispersion), beta-binomial
    (logit-mean with concentration).
  - Missing x handling: latent draws for continuous; marginalized
    discrete x via `log_sum_exp`; shared inference for outcome and
    imputation parameters.
  - Computation: Stan via `instantiate`; vectorized main-model
    likelihood; precomputed id indices; diagnostics (R-hat, ESS,
    divergences, E-BFMI).
- Frequentist comparator (ACML): ODS-only correction using sampling
  probabilities; robust SEs; scope and assumptions.

## 4. Simulation Study Plan

- Objectives: assess bias/RMSE, 95% coverage, interval width, power for
  (*x) and (*{tx}), imputation accuracy for x, cost-adjusted efficiency,
  and computational performance.
- Data generation factors (use
  [`generate_data()`](https://maxdrohde.github.io/bayes2stage/reference/generate_data.md)
  / `generate_data_schild_2019_supp`):
  - N: 500, 1000, 2000; time points M: 3–8; random-effect SDs
    (small/medium/large) and correlations (0, 0.3, 0.6); residual SD.
  - Effects: (*x) and (*{tx}) at 0 (type I), moderate, strong; include
    null for error control.
  - Covariate types: normal; Bernoulli (p=0.2, 0.5); neg-bin (vary
    mean/dispersion); beta-binomial (vary n_trials, concentration).
  - Auxiliary strength: vary gamma coefficients linking z to x
    (weak/medium/strong).
- Design factors:
  - Sampling fractions: 5–40% of subjects; cutoff_low/high (e.g.,
    0.2/0.8); stratum proportions (symmetric vs tail-heavy).
  - Designs: SRS, ODS-intercept, ODS-slope, BDS-intercept, BDS-slope.
  - Misspecification: design on slope when interaction absent/present;
    imputation family misspecified (e.g., fit normal to count x).
- Estimators compared:
  - Bayesian joint model (correct vs misspecified imputation family).
  - ACML (ODS only).
  - Naive complete-case (ignores design) and full-data oracle benchmark.
- Metrics & reporting:
  - Bias, RMSE, interval width, coverage; power/type I error; relative
    efficiency vs oracle.
  - Imputation quality: MAE of x, coverage of x posteriors, calibration
    plots.
  - Cost-adjusted efficiency: precision per unit cost (assign cost
    multiplier to stage-2 x measurement).
  - Computation: wall-clock, divergences, ESS/sec.
- Replication: ≥500 replicates per cell; Cantor seeds for
  reproducibility; summarize via means/quantiles with uncertainty bands.

## 5. Results Presentation

- Design visuals: histograms/density of OLS/BLUP intercepts/slopes with
  cutoff lines; selection heatmaps over time.
- Simulation outputs:
  - Tables of bias/coverage/efficiency by design and x distribution.
  - Power and relative-efficiency curves vs sampling fraction and effect
    size.
  - Forest/ridge plots of (*x) and (*{tx}) posteriors across designs;
    overlay ACML point/intervals.
  - Imputation diagnostics: observed vs imputed x distributions;
    calibration for discrete x; scatter of true vs posterior mean x for
    missing cases.
  - Cost vs precision plots; computational diagnostics (divergence
    rates, ESS/sec).
- Case study: stratum construction plot, fitted trajectories, forest
  plot of key effects, posterior predictive checks for y and x.

## 6. Real-Data Illustration (if available)

- Context and variables; stage-2 design setup (cutoffs, achieved
  sampling proportions).
- Fit Bayesian model matching x type; ACML for ODS where applicable.
- Report effects, uncertainty, and practical interpretation; posterior
  predictive checks for y and x.

## 7. Discussion

- Key findings: when ODS/BDS outperform SRS; value of auxiliary z;
  robustness to imputation misspecification.
- Guidance: choosing cutoffs/strata, recommended sampling fractions vs
  RE variance/M, minimum M for stable BLUPs.
- Limitations: MAR given design variables; linear trajectories;
  baseline-only x; no x measurement error; computational demands at
  larger N.
- Future work: time-varying expensive covariates; nonlinear
  trajectories; informative dropout; adaptive designs; alternative
  priors.

## 8. Reproducibility

- Package: `bayes2stage` (Stan code in `src/stan/`), functions for
  design/fitting
  ([`fit_stan_model()`](https://maxdrohde.github.io/bayes2stage/reference/fit_stan_model.md),
  [`fit_acml_ods()`](https://maxdrohde.github.io/bayes2stage/reference/fit_acml_ods.md)).
- Provide simulation scripts, seeds, session info, and repository;
  suggest container/renv for exact replication.

## 9. Appendices

- Full Stan model listings and priors; ACML derivation; additional
  simulation tables/diagnostics; sensitivity analyses (adapt_delta,
  chains, priors).
