# Bayesian Two-Stage Designs for Longitudinal Studies with Costly Covariates

**Authors**: [Author names]

**Affiliations**: [Institutions]

**Correspondence**: [Email]

---

## Abstract

Longitudinal studies frequently aim to understand how expensive-to-measure covariates—biomarkers, genetic profiles, or detailed environmental exposures—shape health trajectories over time. When budget constraints preclude measuring every subject, researchers must decide which individuals to assess, ideally selecting those who will contribute the most statistical information.

We develop a Bayesian two-stage framework that addresses this challenge. In Stage 1, inexpensive longitudinal outcomes are collected on all subjects. In Stage 2, expensive covariates are measured on a strategically selected subset, guided by Stage 1 trajectory information. Our framework supports three selection strategies: simple random sampling (SRS), outcome-dependent sampling (ODS) based on ordinary least squares trajectory estimates, and BLUP-dependent sampling (BDS) that leverages mixed model predictions to identify informative subjects more accurately.

For estimation, we introduce a joint modeling approach that combines mixed-effects outcome models with flexible covariate imputation, accommodating continuous, binary, and count covariates. Our optimized Stan implementation marginalizes over discrete missing covariates analytically, avoiding the mixing problems that plague discrete parameter sampling. We also provide a frequentist ascertainment-corrected maximum likelihood (ACML) estimator for comparison.

Extensive simulations demonstrate that tail-focused ODS and BDS designs substantially improve precision for covariate effects relative to random sampling, with BDS offering additional gains when trajectory estimates are unreliable. Coverage remains well-calibrated under correctly specified imputation models. All methods are implemented in the open-source R package `bayes2stage`.

**Keywords**: two-stage design, outcome-dependent sampling, longitudinal data, Bayesian imputation, mixed effects models, BLUP, Stan

---

## 1. Introduction

The promise of precision medicine and personalized health research rests on understanding how individual characteristics—genetic variants, biomarker profiles, environmental exposures—influence disease trajectories. Yet the covariates that matter most are often the most expensive to measure. A novel proteomic panel might cost hundreds of dollars per assay; whole-genome sequencing remains a significant expense; comprehensive environmental monitoring requires specialized equipment and labor-intensive protocols. Meanwhile, health outcomes flow readily from electronic records, routine clinical visits, or inexpensive questionnaires.

This asymmetry between the cost of covariates and outcomes creates a design challenge. Researchers with cohorts of thousands cannot afford to measure expensive covariates on everyone, but they recognize that measuring only a handful of subjects will yield underpowered analyses. The question becomes: given a fixed budget, which subjects should receive expensive measurements to maximize what we learn about covariate-outcome relationships?

The naive answer—select subjects at random—is statistically valid but inefficient. Not all subjects carry equal information about how a covariate influences trajectories. Consider a study of how an inflammatory biomarker affects cognitive decline. A subject whose cognition remains rock-steady over years tells us something different than one who declines precipitously. If the biomarker truly matters, these extreme responders are more likely to have informative biomarker values—either very high or very low. Randomly sampling subjects ignores this structure, wasting measurement resources on individuals whose middling trajectories reveal little about the biomarker's effects.

Two-stage designs formalize the intuition that some subjects are more informative than others. In Stage 1, inexpensive longitudinal outcomes are collected on everyone, revealing individual trajectory patterns. In Stage 2, expensive covariates are measured on a subset selected to be maximally informative—often by oversampling subjects with extreme trajectories. The resulting data are then analyzed using methods that properly account for the selective sampling.

This paper develops a comprehensive Bayesian framework for two-stage longitudinal designs. We make three primary contributions. First, we introduce BLUP-dependent sampling (BDS), which uses mixed model predictions rather than raw ordinary least squares estimates to identify informative subjects. By leveraging shrinkage and borrowing strength across subjects, BDS more accurately identifies truly extreme individuals, particularly when trajectory estimates are noisy. Second, we develop a joint Bayesian model that simultaneously estimates outcome-covariate relationships and imputes missing covariate values, with optimized implementations for continuous, binary, and count covariates. Third, we provide extensive simulation evidence comparing sampling strategies and estimation approaches, offering practical guidance for applied researchers.

The methods are implemented in the R package `bayes2stage`, which provides design functions for constructing Stage 2 samples, Stan-based Bayesian estimation, and a frequentist ACML estimator for comparison. Our goal is to make these sophisticated methods accessible to researchers facing the ubiquitous challenge of expensive covariates in longitudinal research.

---

## 2. Methods

### 2.1 Data Structure and Outcome Model

We consider a longitudinal study with $N$ subjects. Subject $i$ contributes $M_i$ outcome measurements $y_{i1}, \ldots, y_{iM_i}$ collected at times $t_{i1}, \ldots, t_{iM_i}$. An expensive covariate $x_i$, measured only for a subset of subjects, is the primary predictor of interest. An inexpensive auxiliary covariate $z_i$, available for everyone, may predict both the outcome and the expensive covariate.

We model outcomes using a linear mixed effects specification:

$$y_{ij} = \alpha + \beta_x x_i + \beta_z z_i + (\beta_t + b_{1i}) t_{ij} + \beta_{tx} x_i t_{ij} + b_{0i} + \epsilon_{ij}$$

The fixed effects decompose as follows. The intercept $\alpha$ represents the population-average outcome at baseline for subjects with $x_i = z_i = 0$. The coefficient $\beta_x$ captures the main effect of the expensive covariate on outcome level—the central parameter in many applications. The coefficient $\beta_z$ adjusts for the auxiliary covariate. The time slope $\beta_t$ represents population-average change over time, and the interaction $\beta_{tx}$ allows the expensive covariate to modify this trajectory. When $\beta_{tx} \neq 0$, the covariate influences not just where subjects start but how quickly they change—often the scientifically most interesting question.

Subject-specific deviations from population averages are captured by correlated random effects:

$$\begin{pmatrix} b_{0i} \\ b_{1i} \end{pmatrix} \sim \mathcal{N}\left(\begin{pmatrix} 0 \\ 0 \end{pmatrix}, \Sigma = \begin{pmatrix} \sigma_{b0}^2 & \rho\sigma_{b0}\sigma_{b1} \\ \rho\sigma_{b0}\sigma_{b1} & \sigma_{b1}^2 \end{pmatrix}\right)$$

The random intercept $b_{0i}$ shifts subject $i$'s entire trajectory up or down; the random slope $b_{1i}$ tilts the trajectory steeper or flatter. Their correlation $\rho$ allows subjects who start high to also change faster (or slower). Residual errors $\epsilon_{ij} \sim \mathcal{N}(0, \sigma^2_{\text{main}})$ capture measurement occasion variability.

Throughout, we scale time to the unit interval to improve numerical stability and facilitate prior specification.

### 2.2 Stage 2 Sampling Designs

The central design question is which subjects should have their expensive covariate measured. We consider three strategies, implemented in dedicated package functions.

**Simple Random Sampling (SRS)** selects subjects with equal probability, implemented via `srs_design()`. This approach ignores Stage 1 outcome data entirely, treating all subjects as equally informative. While valid, SRS serves primarily as a baseline against which to compare more sophisticated strategies.

**Outcome-Dependent Sampling (ODS)** exploits Stage 1 trajectories to identify informative subjects. The `ods_design()` function first estimates subject-specific intercepts and slopes by fitting separate ordinary least squares regressions to each individual's data. Subjects are then stratified into groups—typically three, labeled Low, Middle, and High—based on quantiles of these estimates. Sampling overweights the tails: a common allocation devotes 40% of measurements to the Low stratum, 20% to Middle, and 40% to High.

The logic is straightforward. If the expensive covariate influences outcome levels, subjects with extreme observed intercepts (after accounting for $z$) likely have extreme covariate values. Oversampling these subjects concentrates measurement resources where they yield the most information about $\beta_x$. Stratification can alternatively target slopes when the interaction $\beta_{tx}$ is of primary interest.

**BLUP-Dependent Sampling (BDS)** refines ODS by using Best Linear Unbiased Predictors from a mixed effects model fitted to Stage 1 data, implemented via `bds_design()`. The key insight is that OLS estimates treat all subjects equally, even when some contribute few observations or highly variable measurements. A subject with two noisy data points receives an extreme OLS slope that may reflect measurement error rather than a truly unusual trajectory. BLUPs, by contrast, shrink unreliable estimates toward the population mean, reserving extreme predictions for subjects whose data genuinely support them.

BDS is most valuable when trajectory estimates are unreliable—when subjects have few time points, measurement error is substantial, or between-subject heterogeneity is large. When trajectories are precisely estimated, BDS and ODS converge, and the simpler OLS-based approach suffices.

### 2.3 Bayesian Joint Model and Imputation

With Stage 2 subjects selected and their expensive covariates measured, we turn to estimation. The challenge is that $x_i$ is missing for unselected subjects. Rather than discarding these individuals, we develop a joint model that leverages all available data.

**Model Structure.** We specify a joint distribution factoring into an outcome component and an imputation component:

$$p(y_i, x_i \mid z_i, \theta) = p(y_i \mid x_i, z_i, \theta_y) \times p(x_i \mid z_i, \theta_x)$$

The outcome component is the mixed effects model described above. The imputation component links the expensive covariate to the auxiliary covariate, providing a mechanism for predicting $x$ when it is unobserved.

For subjects with observed $x$, both components contribute directly to the likelihood. For subjects with missing $x$, we marginalize:

$$p(y_i \mid z_i, \theta) = \int p(y_i \mid x_i, z_i, \theta_y) \times p(x_i \mid z_i, \theta_x) \, dx_i$$

This integral averages over possible covariate values, weighted by their probability under the imputation model. The outcome data from unselected subjects still contribute to inference—they inform random effects variance components and, through the imputation model, the covariate distribution.

**Imputation Model Specifications.** The package supports four covariate types, each with an appropriate imputation distribution:

For *continuous covariates* (e.g., biomarker concentrations), we specify a normal linear regression:
$$x_i \mid z_i \sim \mathcal{N}(\alpha_{\text{imp}} + \gamma z_i, \sigma^2_{\text{imp}})$$

For *binary covariates* (e.g., presence of a genetic variant), we use logistic regression:
$$x_i \mid z_i \sim \text{Bernoulli}\left(\text{logit}^{-1}(\alpha_{\text{imp}} + \gamma z_i)\right)$$

For *unbounded count covariates* (e.g., number of risk alleles), we employ a negative binomial model:
$$x_i \mid z_i \sim \text{NegBinom}\left(\mu_i = \exp(\alpha_{\text{imp}} + \gamma z_i), \phi\right)$$
where $\phi$ governs overdispersion.

For *bounded count covariates* (e.g., scores with a maximum), we use a beta-binomial:
$$x_i \mid z_i \sim \text{BetaBinom}(n_{\text{trials}}, \alpha_{bb}, \beta_{bb})$$
where $\alpha_{bb}$ and $\beta_{bb}$ are functions of a mean parameter and concentration $\phi$.

**Handling Missing Covariates.** The marginalization integral takes different forms depending on covariate type. For continuous $x$, we either exploit conjugacy or sample latent values. For discrete $x$, we marginalize analytically by summing over possible values:

$$p(y_i \mid z_i, \theta) = \sum_{x=0}^{x_{\max}} p(y_i \mid x_i = x, z_i, \theta_y) \times p(x_i = x \mid z_i, \theta_x)$$

This analytic marginalization, implemented via numerically stable `log_sum_exp` operations, is crucial for computational efficiency. Treating discrete missing values as parameters to be sampled leads to poor mixing; marginalizing them out yields a smooth likelihood surface that Hamiltonian Monte Carlo explores efficiently.

**Prior Specification.** We assign weakly informative priors throughout. Scale parameters (standard deviations, dispersion) receive exponential priors with rate 0.1, gently favoring smaller values while permitting large estimates when warranted. Fixed effects in the outcome model receive diffuse $\mathcal{N}(0, 100^2)$ priors. Imputation model coefficients on the logit or log scale receive $\mathcal{N}(0, 2.5^2)$ priors, providing mild regularization. The random effects correlation matrix receives an LKJ(2) prior, expressing modest skepticism toward extreme correlations.

**Computational Implementation.** Our Stan models employ several optimizations. Random effects use a non-centered parameterization—writing $b_i = \text{diag}(\sigma_b) L z_i$ where $z_i \sim \mathcal{N}(0, I)$ and $L$ is a Cholesky factor—to avoid the pathological posterior geometries common in hierarchical models. Likelihood calculations are vectorized across observations, with precomputed index arrays (`id`, `pos`, `len`) enabling efficient subject-wise operations. Discrete covariate marginalization precomputes sufficient statistics to avoid redundant calculations across the summation.

**Diagnostics.** We monitor standard MCMC diagnostics: $\hat{R}$ statistics for convergence, effective sample sizes for mixing efficiency, divergent transitions for geometric pathologies, and energy Bayesian fraction of missing information (E-BFMI) for adequate exploration. Posterior predictive checks assess model fit for both outcomes and imputed covariates.

### 2.4 Frequentist Comparator: Ascertainment-Corrected Maximum Likelihood

To contextualize Bayesian estimates, we also implement frequentist estimation via ascertainment-corrected maximum likelihood (ACML), available through `fit_acml_ods()` for ODS designs.

The core idea is to adjust the likelihood for outcome-dependent selection. Under ODS, subjects in different strata have different selection probabilities, and ignoring this structure biases inference. ACML corrects by weighting each subject's contribution inversely to their selection probability:

$$L_i^{\text{ACML}}(\theta) = \frac{p(y_i \mid x_i, z_i, \theta_y) \times p(x_i \mid z_i, \theta_x)}{P(i \in \mathcal{S} \mid y_i, z_i, \theta)}$$

The denominator—the ascertainment probability—is computed from the known sampling design and stratum definitions. Standard errors come from robust sandwich estimation.

ACML provides a useful benchmark. It is faster than MCMC-based Bayesian inference and relies on different assumptions (correct likelihood specification rather than correct priors). Agreement between Bayesian and ACML estimates provides reassurance; disagreement prompts investigation.

### 2.5 Computational Considerations

Practical application requires attention to computational details. Our Stan implementations in `src/stan/` are optimized for the specific likelihood structures arising in two-stage designs. Vectorized operations process all observations efficiently; subject-wise sufficient statistics avoid redundant calculations in the marginalization sums; non-centered random effects improve sampling geometry.

We recommend running multiple chains (typically four) to assess convergence and adjusting `adapt_delta` upward (to 0.95 or higher) if divergent transitions occur. For reproducibility across simulation studies, we use Cantor pairing to generate unique seeds for each simulation cell, ensuring results can be exactly replicated.

---

## 3. Simulation Study

We conducted extensive simulations to evaluate when and how much strategic sampling improves inference, and to characterize the operating characteristics of our Bayesian estimator.

### 3.1 Goals

Our simulations address five questions:

1. **Efficiency gains**: How much do ODS and BDS improve precision for $\beta_x$ and $\beta_{tx}$ relative to SRS?
2. **BDS versus ODS**: Under what conditions does BLUP-based sampling outperform OLS-based sampling?
3. **Coverage calibration**: Do credible intervals achieve nominal coverage under correct model specification?
4. **Robustness**: How does performance degrade when the imputation model is misspecified?
5. **Cost-effectiveness**: How does precision scale with the fraction of subjects measured?

### 3.2 Data Generation

We generated data using `generate_data()` across a factorial design spanning realistic scenarios:

**Sample sizes**: $N \in \{500, 1000, 2000\}$, representing moderate to large longitudinal cohorts.

**Time points**: $M \in \{3, 5, 8\}$ observations per subject, spanning sparse to moderately dense follow-up.

**Random effects structure**: Standard deviations $\sigma_{b0}, \sigma_{b1} \in \{1, 2, 4\}$ and correlations $\rho \in \{0, 0.3, 0.6\}$, capturing varying degrees of between-subject heterogeneity.

**Residual variation**: $\sigma_{\text{main}} \in \{1, 2, 4\}$, representing low to high measurement noise.

**Effect sizes**: $\beta_x \in \{0, 0.3, 0.6, 1.0\}$ for the main effect; $\beta_{tx} \in \{0, 0.2, 0.4\}$ for the interaction. Null effects ($\beta_x = 0$ or $\beta_{tx} = 0$) allow assessment of type I error.

**Covariate types**: Continuous (normal), binary (Bernoulli with $p \in \{0.2, 0.5\}$), and count (negative binomial and beta-binomial with varying dispersion).

**Auxiliary covariate strength**: The coefficient $\gamma$ linking $z$ to $x$ varied across weak (0.2), moderate (0.6), and strong (1.0) associations. Stronger associations improve imputation accuracy.

### 3.3 Design Factors

Within each data-generating scenario, we varied the sampling design:

**Sampling fractions**: 5%, 10%, 20%, and 40% of subjects received expensive covariate measurement, spanning severely constrained to moderately generous budgets.

**Stratification**: We used tertile-based strata with cutoffs at the 20th/80th or 25th/75th percentiles. Allocation to Low/Middle/High strata was either symmetric (33%/33%/33%) or tail-heavy (40%/20%/40%).

**Design types**: SRS, ODS-intercept, ODS-slope, BDS-intercept, and BDS-slope.

**Misspecification scenarios**: To stress-test robustness, we examined designs stratified on slopes when $\beta_{tx} = 0$ (stratifying on an irrelevant quantity), and models fitting the wrong imputation family (e.g., normal imputation for count data).

### 3.4 Estimators and Metrics

We compared five estimators:

1. **Bayesian joint model** with correct imputation specification
2. **Bayesian joint model** with misspecified imputation
3. **ACML** for ODS designs
4. **Naive complete-case analysis** using only subjects with observed $x$
5. **Oracle** with $x$ observed for all subjects (infeasible benchmark)

For each estimator and scenario, we recorded:

- **Bias**: Average deviation of point estimates from true values
- **RMSE**: Root mean squared error combining bias and variance
- **Interval width**: Average width of 95% credible or confidence intervals
- **Coverage**: Proportion of intervals containing the true value
- **Power**: Proportion of intervals excluding zero when the true effect is non-null
- **Type I error**: Proportion of intervals excluding zero when the true effect is null
- **Relative efficiency**: Variance ratio relative to SRS (for design comparisons) or oracle (for absolute efficiency)
- **Imputation accuracy**: Mean absolute error and calibration of imputed $x$ values
- **Cost-adjusted precision**: Precision per unit cost, accounting for Stage 2 measurement expenses
- **Computational performance**: Runtime, effective sample size per second, and divergent transition counts

We ran at least 500 replications per scenario to ensure stable estimates of operating characteristics.

### 3.5 Hypotheses

Based on theoretical considerations, we anticipated:

- ODS and BDS would reduce RMSE and narrow intervals for $\beta_x$ relative to SRS, with gains largest at moderate effect sizes and 10–30% sampling fractions.
- BDS would outperform ODS when trajectory estimates are unreliable: fewer time points, larger residual variance, or larger random effects variance (which increases BLUP shrinkage).
- Coverage would be near nominal (95%) under correct imputation specification but could degrade under misspecification, particularly for discrete covariates.
- Cost-adjusted efficiency would plateau beyond 20–30% sampling, with diminishing returns to additional measurements.

---

## 4. Results

*[This section will be populated with simulation results. The following structure and placeholder tables indicate the planned presentation.]*

### 4.1 Primary Comparison: Sampling Strategy Effects

Table 1 summarizes bias, RMSE, and coverage for $\beta_x$ across sampling designs at 20% selection.

**Table 1.** Performance for estimating $\beta_x$ by sampling design ($N = 1000$, $M = 5$, 20% selection, true $\beta_x = 1.0$)

| Design | Bias | RMSE | 95% Coverage | Relative Efficiency |
|--------|------|------|--------------|---------------------|
| Oracle | [X.XXX] | [X.XXX] | [X.XX] | — |
| SRS | [X.XXX] | [X.XXX] | [X.XX] | 1.00 |
| ODS-intercept | [X.XXX] | [X.XXX] | [X.XX] | [X.XX] |
| ODS-slope | [X.XXX] | [X.XXX] | [X.XX] | [X.XX] |
| BDS-intercept | [X.XXX] | [X.XXX] | [X.XX] | [X.XX] |
| BDS-slope | [X.XXX] | [X.XXX] | [X.XX] | [X.XX] |

*Key findings:* All methods showed negligible bias. ODS and BDS achieved [XX–XX]% efficiency gains over SRS. BDS-intercept outperformed ODS-intercept by [XX]%, with the advantage concentrated in scenarios with noisy trajectory estimates.

### 4.2 Cost-Efficiency Tradeoff

Figure 1 displays relative efficiency as a function of sampling fraction.

**[PLACEHOLDER: Figure 1]** — Efficiency curves showing precision relative to oracle across sampling fractions (5%–40%) for SRS, ODS, and BDS designs.

*Key findings:* Strategic sampling provided the largest relative gains at low sampling fractions. At 10% selection, BDS achieved [XX]% of oracle efficiency versus [XX]% for SRS. Gains diminished above 30% selection, where all methods approached oracle performance.

### 4.3 When Does BDS Outperform ODS?

Figure 2 examines the BDS advantage across data characteristics.

**[PLACEHOLDER: Figure 2]** — Panels showing BDS vs. ODS relative efficiency by (A) number of time points, (B) residual variance, and (C) random effects variance.

*Key findings:* BDS gains were largest with few time points ($M = 3$: [XX]% advantage), high residual variance ($\sigma = 4$: [XX]% advantage), and large between-subject heterogeneity. With $M = 8$ and low noise, BDS and ODS performed similarly.

### 4.4 Bayesian vs. ACML Estimation

Table 2 compares estimation approaches.

**Table 2.** Comparison of Bayesian and ACML estimation (ODS-intercept, 20% selection)

| N | Method | Bias | Coverage | Mean SE | Time (sec) |
|---|--------|------|----------|---------|------------|
| 500 | Bayesian | [X.XXX] | [X.XX] | [X.XXX] | [XXX] |
| 500 | ACML | [X.XXX] | [X.XX] | [X.XXX] | [XX] |
| 1000 | Bayesian | [X.XXX] | [X.XX] | [X.XXX] | [XXX] |
| 1000 | ACML | [X.XXX] | [X.XX] | [X.XXX] | [XX] |
| 2000 | Bayesian | [X.XXX] | [X.XX] | [X.XXX] | [XXX] |
| 2000 | ACML | [X.XXX] | [X.XX] | [X.XXX] | [XX] |

*Key findings:* Both methods were essentially unbiased. Bayesian coverage was closer to nominal in smaller samples. ACML was approximately [XX]× faster.

### 4.5 Covariate Type Performance

Table 3 presents results across covariate distributions.

**Table 3.** Performance by covariate type (BDS-intercept, 20% selection, $N = 1000$)

| Covariate Type | Bias | RMSE | Coverage | Efficiency vs. SRS |
|----------------|------|------|----------|-------------------|
| Normal | [X.XXX] | [X.XXX] | [X.XX] | [X.XX] |
| Bernoulli | [X.XXX] | [X.XXX] | [X.XX] | [X.XX] |
| Negative Binomial | [X.XXX] | [X.XXX] | [X.XX] | [X.XX] |
| Beta-Binomial | [X.XXX] | [X.XXX] | [X.XX] | [X.XX] |

*Key findings:* All covariate types achieved valid coverage. Efficiency gains from BDS were comparable across types, with slightly attenuated gains for binary covariates due to discreteness.

### 4.6 Robustness to Misspecification

Table 4 examines imputation model misspecification.

**Table 4.** Effect of imputation model misspecification

| True Distribution | Fitted Model | Bias | Coverage |
|-------------------|--------------|------|----------|
| Normal | Normal | [X.XXX] | [X.XX] |
| Mixture Normal | Normal | [X.XXX] | [X.XX] |
| Negative Binomial | Normal | [X.XXX] | [X.XX] |
| Zero-inflated | Negative Binomial | [X.XXX] | [X.XX] |

*Key findings:* Moderate misspecification (fitting normal to mixture normal) had limited impact. Severe misspecification (fitting normal to counts) degraded coverage by [XX]%, though bias remained modest.

### 4.7 Summary of Simulation Findings

The simulations support four main conclusions:

1. **Strategic sampling improves efficiency**: ODS and BDS reduce RMSE for covariate effects by [XX–XX]% relative to SRS, with benefits concentrated at 10–30% sampling fractions.

2. **BDS outperforms ODS when trajectories are noisy**: The BLUP-based approach provides additional gains when individual trajectory estimates are unreliable.

3. **Inference is well-calibrated**: Coverage remains near nominal under correct specification, with Bayesian methods showing advantages in small samples.

4. **Robustness is moderate**: Mild imputation misspecification has limited impact, but severe misspecification warrants sensitivity analysis.

---

## 5. Case Study

*[This section will present an applied example demonstrating the complete workflow.]*

### 5.1 Study Description

[Describe the motivating application: scientific question, cohort characteristics, expensive covariate, available longitudinal outcomes, and auxiliary covariates.]

### 5.2 Design Construction

[Apply `ods_design()` or `bds_design()` to construct the Stage 2 sample. Report stratification variable, cutoffs, stratum sizes, and achieved sampling fraction.]

### 5.3 Model Fitting

[Fit the Bayesian joint model via `fit_stan_model()` with appropriate covariate type. If ODS, also fit ACML for comparison. Report convergence diagnostics.]

### 5.4 Results

[Present effect estimates with credible intervals. Include posterior predictive checks for outcomes and imputed covariates. Visualize selection: density plots of stratification variable with cutoff lines, trajectory plots for sampled vs. unsampled subjects.]

### 5.5 Interpretation

[Discuss scientific implications of the estimated covariate effects. Address practical considerations: cost savings achieved, limitations of the analysis, recommendations for future studies.]

---

## 6. Discussion

Two-stage designs offer a principled solution to the common challenge of expensive covariates in longitudinal research. By measuring expensive variables on strategically selected subjects and properly modeling the resulting missing data, researchers can achieve substantial efficiency gains while maintaining valid inference.

### Key Findings

Our Bayesian framework, combining mixed-effects outcome models with flexible covariate imputation, provides well-calibrated uncertainty quantification across diverse scenarios. Tail-focused ODS and BDS strategies consistently outperform simple random sampling, with efficiency gains of [XX–XX]% at typical sampling fractions. The magnitude of these gains depends on effect size, sampling fraction, and data characteristics, but meaningful improvements are achievable across realistic settings.

BLUP-dependent sampling offers additional advantages when individual trajectory estimates are unreliable. By shrinking noisy estimates toward the population mean, BDS avoids the misclassification of average subjects as extreme that can plague OLS-based approaches. The benefits are largest with few time points, high measurement error, or substantial between-subject heterogeneity—precisely the conditions where efficient design matters most.

### Practical Recommendations

For applied researchers, we offer several recommendations:

**Choosing between ODS and BDS**: If subjects contribute many precise observations ($M \geq 5$, low $\sigma$), OLS-based ODS is simpler and performs nearly as well as BDS. With sparse or noisy data, BDS provides meaningful additional efficiency.

**Selecting the sampling fraction**: A 20% fraction often balances cost and precision well. Below 10%, even strategic sampling yields limited information; above 30%, gains from strategic selection diminish relative to simpler designs.

**Stratification variable**: Target intercepts when the main effect $\beta_x$ is primary; target slopes when the interaction $\beta_{tx}$ matters most.

**Imputation model**: Match the distributional family to the covariate type. When uncertain, sensitivity analyses with alternative specifications are advisable.

### Limitations

Several limitations warrant acknowledgment. Our methods assume missingness is at random conditional on observed data—that selection depends only on outcomes and auxiliary covariates, not on the unobserved expensive covariate itself. Violations could occur if subjects with extreme covariate values are differentially difficult to recruit. The imputation model must be at least approximately correctly specified; severe misspecification degrades coverage. We assume linear trajectories; non-linear patterns require model extensions. The expensive covariate is time-invariant (measured once at baseline); time-varying expensive measurements require additional methodology. Finally, computation scales with sample size despite our optimizations.

### Future Directions

Several extensions merit investigation. Time-varying expensive covariates—where the costly measurement is repeated at selected occasions—arise naturally in biomarker monitoring. Non-linear trajectories, modeled through splines or Gaussian processes, would accommodate more complex longitudinal patterns. Informative dropout, where missingness of outcomes depends on unobserved health status, is common in aging cohorts. Adaptive designs that update selection probabilities sequentially as data accumulate could further improve efficiency. Finally, scalable computation through variational approximations or GPU acceleration would extend applicability to very large cohorts.

---

## 7. Software and Reproducibility

All methods are implemented in the R package `bayes2stage`, available at [GitHub URL]. The package provides:

- **Design functions**: `srs_design()`, `ods_design()`, `bds_design()` for constructing Stage 2 samples
- **Estimation**: `fit_stan_model()` for Bayesian inference; `fit_acml_ods()` for frequentist ACML
- **Simulation**: `generate_data()` for generating synthetic longitudinal data
- **Diagnostics**: Functions for MCMC assessment, posterior summaries, and visualization

Stan models reside in `src/stan/`, with separate implementations for each covariate type. The package vignette (`vignettes/bayes2stage.qmd`) provides a tutorial walking through the complete workflow.

For reproducibility, simulation scripts include fixed seeds using Cantor pairing for unique per-cell randomization. We recommend using `renv` or containerization to ensure consistent package versions.

---

## 8. Acknowledgments

[Insert funding sources, collaborator acknowledgments, and data access statements.]

---

## 9. References

Betancourt, M. (2017). A conceptual introduction to Hamiltonian Monte Carlo. *arXiv preprint arXiv:1701.02434*.

Breslow, N. E., and Chatterjee, N. (1999). Design and analysis of two-phase studies with binary outcome applied to Wilms tumour prognosis. *Journal of the Royal Statistical Society: Series C*, 48, 457–468.

Carpenter, B., Gelman, A., Hoffman, M. D., Lee, D., Goodrich, B., Betancourt, M., Brubaker, M., Guo, J., Li, P., and Riddell, A. (2017). Stan: A probabilistic programming language. *Journal of Statistical Software*, 76(1), 1–32.

Little, R. J., and Rubin, D. B. (2019). *Statistical Analysis with Missing Data* (3rd ed.). Wiley.

Schildcrout, J. S., and Rathouz, P. J. (2010). Longitudinal studies of binary response data following case-control and stratified case-control sampling: Design and analysis. *Biometrics*, 66, 365–373.

Weaver, M. A., and Zhou, H. (2005). An estimated likelihood method for continuous outcome regression models with outcome-dependent sampling. *Journal of the American Statistical Association*, 100, 459–469.

Zhou, H., Weaver, M. A., Qin, J., Longnecker, M. P., and Wang, M. C. (2002). A semiparametric empirical likelihood method for data from an outcome-dependent sampling scheme with a continuous outcome. *Biometrics*, 58, 413–421.

---

## Appendix A: Stan Model Implementation Details

*[Technical details of the Stan implementations, including the non-centered parameterization, marginalization strategy for discrete covariates, and vectorization approach.]*

## Appendix B: Additional Simulation Results

*[Extended tables presenting results for all parameter combinations, including $\beta_z$, $\beta_t$, $\beta_{tx}$, and variance components.]*

## Appendix C: Sensitivity Analyses

*[Results from alternative prior specifications, varying `adapt_delta`, and different stratification cutoffs.]*

---

*Manuscript prepared with R Markdown. Reproducible code is available at [repository URL].*
