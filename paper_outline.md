# Paper Outline: Two-Stage Bayesian Methods for Cost-Effective Longitudinal Studies with Expensive Covariates

## Working Title Options
1. "Bayesian Two-Stage Designs for Longitudinal Studies: Efficient Estimation When Covariates Are Expensive to Measure"
2. "BLUP-Dependent Sampling in Longitudinal Studies: A Bayesian Imputation Approach"
3. "Cost-Effective Longitudinal Study Design Using Outcome-Dependent Sampling and Bayesian Imputation"

---

## 1. Introduction and Motivation

### 1.1 The Problem
- Many longitudinal studies require measurement of expensive covariates (biomarkers, genetic tests, environmental exposures)
- Measuring all subjects is often cost-prohibitive
- **Key question**: How can we design studies that maximize statistical efficiency while minimizing measurement costs?

### 1.2 Motivating Examples
- **Biomarker studies**: Expensive assays (e.g., proteomics, metabolomics) in longitudinal cohorts
- **Genetic studies**: Whole-genome sequencing costs on subset of epidemiological cohort
- **Environmental health**: Personal exposure monitoring (e.g., air pollution sensors) is labor-intensive
- **Clinical trials**: Collecting expensive imaging or invasive secondary endpoints

### 1.3 Existing Approaches
- Simple random sampling (SRS) of subjects for expensive measurement
- Case-control designs (for binary outcomes)
- Outcome-dependent sampling (ODS) for cross-sectional studies
- **Gap**: Limited methods for longitudinal data with continuous outcomes and trajectory-based sampling

### 1.4 Our Contribution
- Extension of ODS to longitudinal settings using trajectory characteristics
- Novel BLUP-dependent sampling (BDS) that leverages mixed model estimates
- Bayesian imputation framework that jointly models outcomes and missing covariates
- Comparison with ascertainment-corrected maximum likelihood (ACML)
- Flexible framework supporting different covariate distributions (normal, binary, count)

---

## 2. Methods

### 2.1 Study Design and Notation

#### 2.1.1 Data Structure
- N subjects, each with $n_i$ longitudinal observations
- $y_{ij}$: outcome for subject $i$ at time $t_{ij}$
- $x_i$: expensive covariate (measured on subset)
- $z_i$: inexpensive covariate (measured on all)

#### 2.1.2 Two-Stage Design
- **Stage 1**: Collect $(y_{ij}, z_i, t_{ij})$ on all N subjects
- **Stage 2**: Select $n_2 < N$ subjects for expensive covariate measurement based on Stage 1 data

### 2.2 Outcome Model

Present the mixed effects model:
$$y_{ij} = \alpha + \beta_x x_i + \beta_z z_i + (\beta_t + b_{1i}) t_{ij} + \beta_{xt} x_i t_{ij} + b_{0i} + \epsilon_{ij}$$

where:
- $(b_{0i}, b_{1i})' \sim N(0, \Sigma_b)$: correlated random effects
- $\epsilon_{ij} \sim N(0, \sigma^2)$: residual error

### 2.3 Sampling Strategies

#### 2.3.1 Simple Random Sampling (SRS)
- Baseline comparator
- Select $n_2$ subjects uniformly at random

#### 2.3.2 Outcome-Dependent Sampling (ODS)
- Estimate subject trajectories using OLS on Stage 1 data
- Compute empirical BLUPs: $\hat{b}_{0i}$ (intercept) or $\hat{b}_{1i}$ (slope)
- Stratify into L (low), M (middle), H (high) tertiles
- Oversample from tails (e.g., 40% L, 20% M, 40% H)

#### 2.3.3 BLUP-Dependent Sampling (BDS)
- Fit mixed effects model to Stage 1 data (ignoring $x$)
- Extract conditional BLUPs accounting for correlation structure
- Stratify and sample as in ODS
- **Hypothesis**: More efficient than ODS due to shrinkage

### 2.4 Estimation Approaches

#### 2.4.1 Bayesian Imputation
- Joint model: $p(y, x | z, \theta) = p(y | x, z, \theta_y) \times p(x | z, \theta_x)$
- Imputation model: $x_i \sim f(x | z_i, \theta_x)$ with distribution-specific forms
- For missing $x_i$: marginalize over plausible values
- Stan implementation with efficient computation strategies

#### 2.4.2 Ascertainment-Corrected Maximum Likelihood (ACML)
- Correct for biased sampling via ascertainment probability
- $L_{ACML} = \prod_{i \in S} \frac{p(y_i | x_i) p(x_i | z_i)}{AC_i}$
- Robust sandwich variance estimation

### 2.5 Imputation Model Specifications
- **Normal**: $x_i \sim N(\alpha_{imp} + \gamma z_i, \sigma_{imp}^2)$
- **Bernoulli**: $x_i \sim Bernoulli(\text{logit}^{-1}(\alpha_{imp} + \gamma z_i))$
- **Negative Binomial**: $x_i \sim NegBinom(\mu = \exp(\alpha_{imp} + \gamma z_i), \phi)$
- **Beta-Binomial**: For bounded counts with overdispersion

---

## 3. Simulation Studies

### 3.1 Simulation Design Overview

| Factor | Levels | Purpose |
|--------|--------|---------|
| Sample size (N) | 500, 1000, 2000 | Assess scalability |
| Selection fraction | 10%, 20%, 30% | Cost-efficiency trade-off |
| Effect size ($\beta_x$) | 0.5, 1.0, 2.0 | Power analysis |
| Covariate distribution | Normal, Binary, Count | Methodological flexibility |
| Sampling strategy | SRS, ODS, BDS | Core comparison |
| Estimation method | Bayesian, ACML | Frequentist vs. Bayesian |

### 3.2 Primary Simulation: Comparing Sampling Strategies

**Objective**: Demonstrate efficiency gains of ODS/BDS over SRS

**Setup**:
- N = 1000 subjects, 5 time points each
- True parameters: $\alpha = 2$, $\beta_x = 1$, $\beta_z = 0.5$, $\beta_t = 0.3$, $\beta_{xt} = 0.2$
- Random effects: $\sigma_{b0} = 1$, $\sigma_{b1} = 0.3$, $\rho = 0.3$
- Select 20% of subjects (n = 200)

**Metrics**:
- Bias: $E[\hat{\beta}_x] - \beta_x$
- RMSE: $\sqrt{E[(\hat{\beta}_x - \beta_x)^2]}$
- Coverage: 95% CI coverage probability
- Relative efficiency: $\text{Var}_{SRS}(\hat{\beta}_x) / \text{Var}_{method}(\hat{\beta}_x)$

**Visualizations**:
- Box plots of parameter estimates across methods
- Efficiency curves as function of selection fraction
- Coverage probability plots

### 3.3 Simulation 2: ODS vs. BDS Comparison

**Objective**: When does BLUP-dependent sampling outperform OLS-based ODS?

**Factors**:
- ICC (intraclass correlation): 0.1, 0.3, 0.5, 0.7
- Number of time points: 3, 5, 10
- Measurement error variance: low, medium, high

**Hypothesis**: BDS advantage increases with:
- Higher ICC (more borrowing across subjects)
- Fewer time points (more shrinkage benefit)
- Higher measurement error (more regularization benefit)

### 3.4 Simulation 3: Bayesian vs. ACML Estimation

**Objective**: Compare frequentist and Bayesian approaches

**Metrics**:
- Bias and efficiency (as above)
- Computational time
- Behavior in small samples (n = 100, 200)
- Sensitivity to prior specification

**Visualization**:
- Side-by-side forest plots of estimates
- Posterior distributions vs. asymptotic confidence intervals

### 3.5 Simulation 4: Covariate Distribution Misspecification

**Objective**: Robustness when imputation model is wrong

**Scenarios**:
- True: Normal → Fitted: Normal (baseline)
- True: Mixture of normals → Fitted: Normal
- True: Log-normal → Fitted: Normal
- True: Zero-inflated Poisson → Fitted: Negative Binomial

**Metrics**: Bias and coverage degradation

### 3.6 Simulation 5: Sensitivity to Sampling Fraction

**Objective**: How low can we go?

**Selection fractions**: 5%, 10%, 15%, 20%, 30%, 50%, 100%

**Key question**: At what point does ODS/BDS efficiency gain outweigh information loss from not measuring all subjects?

### 3.7 Simulation 6: Interaction Effect Focus

**Objective**: Efficiency for detecting covariate-by-time interaction ($\beta_{xt}$)

**Rationale**: Interaction effects often require larger samples; can ODS/BDS help?

**Comparison**: Power curves for detecting non-zero $\beta_{xt}$ across methods

---

## 4. Results Presentation Ideas

### 4.1 Main Results Table

| Method | $\hat{\beta}_x$ Bias | RMSE | 95% Coverage | Rel. Efficiency |
|--------|---------------------|------|--------------|-----------------|
| Full data | — | — | — | 1.00 (reference) |
| SRS | | | | |
| ODS-intercept | | | | |
| ODS-slope | | | | |
| BDS-intercept | | | | |
| BDS-slope | | | | |

### 4.2 Figure Concepts

**Figure 1: Conceptual Diagram**
- Schematic of two-stage design
- Show stratification and sampling from tails
- Intuition for why tail sampling is informative

**Figure 2: Efficiency Gains Plot**
- X-axis: Selection fraction (10%-50%)
- Y-axis: Relative efficiency vs. SRS
- Lines for ODS-int, ODS-slope, BDS-int, BDS-slope
- Horizontal line at 1.0 (SRS reference)

**Figure 3: Bias-Variance Trade-off**
- Panel A: Bias across methods
- Panel B: Variance across methods
- Panel C: RMSE (combined)

**Figure 4: Coverage Probability**
- Nominal 95% coverage line
- Bar plot showing empirical coverage by method

**Figure 5: Posterior Comparison (Bayesian)**
- Single simulation run
- Overlay posterior densities for key parameters
- Compare imputed vs. observed-only analysis

**Figure 6: Computational Comparison**
- Timing benchmarks: Bayesian vs. ACML
- Scaling with N and number of time points

**Figure 7: Robustness Analysis**
- Heatmap: Bias under model misspecification scenarios
- Rows: True distribution; Columns: Fitted distribution

### 4.3 Supplementary Figures

- MCMC diagnostics (trace plots, R-hat values)
- Prior sensitivity analysis
- Additional parameter estimates ($\beta_z$, $\beta_t$, random effects variances)

---

## 5. Discussion Points

### 5.1 Key Findings (anticipated)
- ODS and BDS achieve 20-40% efficiency gains over SRS for main covariate effect
- BDS outperforms ODS when ICC is high or time points are sparse
- Bayesian imputation provides valid inference with appropriate uncertainty quantification
- ACML offers computational advantages for large datasets

### 5.2 Practical Recommendations
- When to use ODS vs. BDS
- Recommended minimum selection fractions
- Choice of stratification variable (intercept vs. slope)
- Handling of covariate distribution uncertainty

### 5.3 Limitations
- Assumes MAR (missing at random) conditional on observed data
- Requires specification of imputation model
- Computational demands of Bayesian approach
- Sensitivity to stratification boundary choices

### 5.4 Extensions and Future Work
- Multivariate expensive covariates
- Time-varying expensive covariates
- Non-linear trajectories
- Survival outcomes
- Adaptive designs (sequential sampling)

---

## 6. Software Section

### 6.1 Package Description
- R package `bayes2stage` available on GitHub
- Depends on Stan (via `instantiate`) for Bayesian computation
- Reproducible examples in vignette

### 6.2 Code Availability
- Simulation code in supplementary materials
- Documented API with examples

---

## 7. Appendix Ideas

### A. Technical Details
- Full Stan model specification
- Derivation of ACML likelihood
- Proof of asymptotic properties (if applicable)

### B. Additional Simulation Results
- Full tables with all parameters
- Sensitivity analyses
- Convergence diagnostics

### C. Real Data Application (Optional)
- If available: Apply to motivating dataset
- Compare sampling strategies retrospectively

---

## Timeline and Priority

### High Priority Simulations
1. Primary comparison (SRS vs. ODS vs. BDS) - **Core result**
2. Selection fraction sensitivity - **Practical guidance**
3. Bayesian vs. ACML comparison - **Method validation**

### Medium Priority
4. ODS vs. BDS detailed comparison
5. Covariate distribution variants

### Lower Priority (if time permits)
6. Model misspecification robustness
7. Interaction effect power analysis

---

## Key Messages to Convey

1. **Cost savings are achievable**: Measuring 20% of subjects with smart sampling can achieve 80%+ efficiency of full data
2. **Trajectory-based sampling is intuitive**: Oversampling extreme responders is informative for covariate effects
3. **Bayesian imputation handles uncertainty properly**: Unlike complete-case analysis, accounts for missing data mechanism
4. **BDS is novel and often superior**: Leveraging mixed model machinery improves on simple OLS-based ODS
5. **Practical software is available**: `bayes2stage` package makes methods accessible

---

## References to Include

- Outcome-dependent sampling: Zhou et al., Biometrics literature
- Two-phase designs: Breslow & Chatterjee
- Bayesian missing data: Little & Rubin
- Stan and HMC: Carpenter et al., Stan Development Team
- Mixed effects models: Laird & Ware, Verbeke & Molenberghs
