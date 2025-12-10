# Bayesian Two-Stage Designs for Longitudinal Studies with Costly Covariates

*Draft for development; replace placeholder items (data source, numbers) with study specifics.*

## Abstract
Two-stage designs are attractive when a key covariate (e.g., biomarker or omics
assay) is costly but longitudinal outcomes are widely available. We propose a
Bayesian joint modeling framework that combines mixed-effects outcome models
with covariate imputation, enabling efficient inference when the costly
covariate is measured on a subsample selected by outcome- or BLUP-dependent
rules. The `bayes2stage` package implements unified Stan models for continuous,
binary, and count covariates (negative binomial, beta-binomial), optimized for
computational efficiency. We compare outcome-dependent sampling (ODS),
BLUP-dependent sampling (BDS), and simple random sampling (SRS) across a broad
set of simulation scenarios, and benchmark against an ascertainment-corrected
maximum likelihood (ACML) estimator. Results show that tail-focused ODS/BDS
strategies yield substantial precision gains for covariate effects at modest
costs, with well-calibrated uncertainty under correct imputation families. We
illustrate the approach in a case study and provide open-source code for full
reproducibility.

**Keywords:** two-stage design; outcome-dependent sampling; BLUP-dependent
sampling; missing covariates; mixed-effects models; Bayesian imputation; Stan;
ascertainment correction.

## 1. Introduction
Measuring expensive biomarkers or genetic assays on every participant is often
infeasible. In longitudinal studies, repeated outcomes and inexpensive
covariates are typically observed on all subjects. Two-stage designs exploit
this structure: stage 1 collects inexpensive measures on everyone; stage 2
selects a subsample for the costly covariate, often oversampling subjects with
extreme trajectories to maximize information.

Existing approaches either ignore stage-2 sampling in downstream inference
(leading to bias or inefficiency) or rely on frequentist ascertainment
corrections specialized to particular designs. Bayesian methods can jointly
model the outcome and missing covariate, propagate uncertainty, and flexibly
handle diverse covariate types, but practical, optimized implementations are
scarce. We develop a unified Bayesian framework for two-stage longitudinal
studies and provide an open-source implementation in `bayes2stage`. The
framework supports SRS, ODS (intercept/slope strata from OLS fits), and BDS
(intercept/slope strata from mixed-model BLUPs). It handles continuous, binary,
and bounded/unbounded count covariates with optimized Stan code that
marginalizes discrete covariates when possible. We benchmark against an
ascertainment-corrected maximum likelihood (ACML) estimator and provide design
guidance through simulations and a case study.

## 2. Methods

### 2.1 Data structure and outcome model
Subjects \(i=1,\dots,N\) are observed at times \(j=1,\dots,M_i\) (often a common
M). The outcome model is
\[
 y_{ij} = \alpha + \beta_x x_i + \beta_z z_i + (\beta_t + b_{1i}) t_{ij}
          + \beta_{tx} x_i t_{ij} + b_{0i} + \epsilon_{ij},
\]
with correlated random intercept/slope \((b_{0i}, b_{1i}) \sim \mathcal{N}(0, \Sigma)
\) and residuals \(\epsilon_{ij} \sim \mathcal{N}(0, \sigma_{\text{main}}^2)\). Time
is scaled 0–1. Covariate z is inexpensive and observed for all subjects; x is
expensive, measured only for the stage-2 sample.

### 2.2 Stage-2 sampling strategies
- **SRS:** random subsample of subjects via `srs_design()`.
- **ODS:** define strata on OLS intercepts or slopes (`ods_design()`), using
  quantiles (cutoff_low/high) and prespecified stratum proportions; oversample
  tails to emphasize informative trajectories.
- **BDS:** define strata on BLUP intercepts or slopes from a mixed-effects model
  (`bds_design()`); requires sufficient time points for stable BLUPs and is
  sensitive to random-effect variance magnitude.

### 2.3 Bayesian joint model with imputation
We jointly model y and x. Implementation uses optimized Stan code in
`src/stan/`.

- **Outcome component:** non-centered random effects with Cholesky factor
  \(L_{re}\) and SDs \(\sigma_{re}\); fixed effects \(\beta_x, \beta_z, \beta_t,
  \beta_{tx}\); Gaussian residual variance \(\sigma_{\text{main}}^2\). Likelihood
  is vectorized; subject indices (`id`, `pos`, `len`) broadcast random effects
  efficiently.
- **Imputation component (subject-level x):**
  - Normal: linear regression with SD \(\sigma_{\text{imp}}\).
  - Bernoulli: logit link.
  - Negative binomial: log mean with dispersion \(\phi\).
  - Beta-binomial: logit mean with concentration \(\phi\) and known trials
    \(n_{\text{trials}}\).
- **Missing x handling:** continuous x uses latent draws; discrete x are
  marginalized with `log_sum_exp` to avoid discrete parameters and improve
  mixing.
- **Priors:** exponential on scale parameters (rate 0.1), weakly informative
  normals on fixed effects (mean 0, sd 100 for main effects; sd 2.5 on
  imputation logits), LKJ prior on correlation via `lkj_corr_cholesky`. These
  choices balance weak informativeness with computational stability.
- **Diagnostics:** R-hat, ESS, divergences, E-BFMI; graphical PPCs for y and x;
  calibration of imputed x.

### 2.4 Frequentist comparator (ACML)
For ODS designs, `fit_acml_ods()` implements ascertainment-corrected ML. Known
sampling probabilities by stratum correct the likelihood; robust standard errors
provide inference. This serves as a comparator for the Bayesian approach under
ODS (intercept/slope) designs.

### 2.5 Computational considerations
- Stan models leverage precomputed subject indices (`pos`, `len`) and
  vectorization to reduce gradient evaluations.
- Discrete x models integrate out x to avoid discrete parameters.
- Parallel chains via `parallel_chains`; adapt_delta tuned to reduce divergences.
- Seeds via Cantor pairing for reproducibility across simulation cells.

## 3. Simulation Study

### 3.1 Objectives
Assess:
1) Bias/RMSE and 95% coverage for \(\beta_x\) and \(\beta_{tx}\);
2) Interval width and power/type I error;
3) Imputation quality for x (MAE, calibration);
4) Cost-adjusted precision;
5) Computational performance (runtime, ESS/sec, divergences).

### 3.2 Data-generating mechanisms
Use `generate_data()` and `generate_data_schild_2019_supp` variants:
- **Sample sizes:** N ∈ {500, 1000, 2000}; time points M ∈ {3, 5, 8}.
- **Random effects:** SDs {1, 2, 4}; correlations {0, 0.3, 0.6}; residual SD
  {1, 2, 4} (scale to scenario).
- **Effects:** \(\beta_x\) ∈ {0, 0.3, 0.6, 1.0}; \(\beta_{tx}\) ∈ {0, 0.2, 0.4}.
- **Covariate types:** normal; Bernoulli (p=0.2, 0.5); negative binomial
  (vary mean/dispersion); beta-binomial (vary trials/concentration).
- **Auxiliary strength:** gamma coefficients linking z→x at weak/medium/strong
  (e.g., gamma1 ∈ {0.2, 0.6, 1.0}).

### 3.3 Design factors
- Sampling fractions: 5%, 10%, 20%, 40% of subjects measured for x.
- Strata cutoffs: (0.2, 0.8) and (0.25, 0.75); stratum proportions symmetric vs
  tail-heavy (e.g., 0.4/0.2/0.4 low/mid/high).
- Designs: SRS, ODS-intercept, ODS-slope, BDS-intercept, BDS-slope.
- Misspecification: design on slope when \(\beta_{tx}=0\); fit wrong imputation
  family (e.g., normal to count x) to test robustness.

### 3.4 Estimators
- Bayesian joint model (correctly specified) for each x family.
- Bayesian joint model with imputation misspecification (robustness stress).
- ACML for ODS.
- Naive complete-case ignoring design.
- Oracle (full x observed) for efficiency benchmark.

### 3.5 Metrics and analysis plan
- Compute bias, RMSE, interval width, and coverage for \(\beta_x, \beta_{tx}\) per
  cell; summarize across replicates (≥500 per cell) with means and quantiles.
- Power/type I error from rejection of \(H_0: \beta_x=0\).
- Relative efficiency: variance ratio vs oracle and vs SRS.
- Imputation: MAE of x among missing subjects; calibration of x posteriors
  (coverage of latent draws for discrete x).
- Cost-adjusted precision: precision per cost unit with stage-2 cost multiplier
  (e.g., 10× stage-1).
- Computational: wall-clock, ESS/sec, divergences; flag problematic settings.

### 3.6 Expected qualitative patterns (to verify)
- Tail-focused ODS/BDS should reduce RMSE and interval width for \(\beta_x\),
  especially when effects moderate/large and sampling 10–30% of subjects.
- BDS may outperform ODS when BLUPs are stable (larger M, larger RE variance);
  ODS more robust at small M.
- Misspecified imputation degrades coverage for discrete x; normal fit to counts
  inflates bias when dispersion large.
- Power gains taper beyond ~20–30% sampling; cost-adjusted efficiency plateaus.

## 4. Results (populate with outputs)

### 4.1 Main simulation tables
- Bias/RMSE/coverage for \(\beta_x, \beta_{tx}\) by design and x type.
- Relative efficiency vs SRS and oracle; power/type I error summaries.

### 4.2 Figures
- Efficiency vs sampling fraction curves (per design, per x type).
- Forest/ridge plots of \(\beta_x, \beta_{tx}\) posteriors across designs;
  overlay ACML point/CI.
- Imputation diagnostics: observed vs imputed x distributions; calibration plots
  for discrete x; scatter of true vs posterior mean x among missing cases.
- Cost vs precision trade-offs; computational diagnostics (ESS/sec, divergences).
- Design visuals: density of OLS/BLUP intercepts/slopes with cutoff lines;
  selection heatmaps over time.

### 4.3 Sensitivity analyses
- Vary adapt_delta/step size to assess robustness of HMC diagnostics.
- Alternative priors (tighter/looser scales) for stability and shrinkage effects.
- Impact of weaker z–x association on imputation performance.

## 5. Case Study (insert application)
- **Setting:** describe population, outcome, expensive covariate, and auxiliary z.
- **Design:** construct ODS/BDS strata (intercept or slope), report cutoffs and
  achieved sampling fractions; compare to hypothetical SRS.
- **Models:** fit Bayesian joint model matching x type; fit ACML for ODS if
  applicable.
- **Outputs:** effect estimates with intervals; posterior predictive checks for y
  and x; selection visualizations; forest plot of key effects; trajectory plots
  for sampled vs unsampled subjects.
- **Interpretation:** practical implications of estimated \(\beta_x, \beta_{tx}\);
  cost implications of sampling strategy.

## 6. Discussion
- **Key findings:** tail-focused two-stage designs with joint Bayesian imputation
  improve efficiency for costly covariate effects; gains largest when effects or
  trajectory heterogeneity are moderate/large and auxiliary z is informative.
- **Design guidance:** pick cutoffs to target informative tails; ensure enough
  time points for BLUP stability before using BDS; sample 10–30% when costs are
  high; align design with estimand (use slope-based strata when interaction is
  primary target).
- **Robustness:** correct imputation family matters for discrete x; auxiliary z
  improves calibration; BDS sensitive to M and RE variance; ODS simpler and more
  stable when M is small.
- **Limitations:** assumes MAR given design variables; linear trajectories;
  baseline (time-invariant) x; no measurement error in x; computation grows with
  N despite vectorization.
- **Future work:** time-varying expensive covariates; nonlinear trajectories;
  informative dropout; adaptive or response-adaptive designs; alternative priors
  and variational approximations for scalability.

## 7. Software and Reproducibility
- Package: `bayes2stage`; Stan code in `src/stan/`; design functions
  (`srs_design()`, `ods_design()`, `bds_design()`); Bayesian fitting via
  `fit_stan_model()`; ACML via `fit_acml_ods()`.
- Reproducibility: publish simulation scripts, seeds, session info; recommend
  `renv` or containers; vignette (`vignettes/bayes2stage.qmd`) shows workflow.

## 8. Acknowledgments
[Insert funding, collaborators, and data acknowledgments.]

## 9. References
[Insert citations for two-stage designs, ODS/BDS literature, ACML, Bayesian
missing data, Stan references, and the applied domain.] 
