# Practical Guide to Two-Stage Designs

## Executive Summary

This guide provides practical recommendations for applied researchers considering two-stage designs when covariates are expensive to measure in longitudinal studies.

---

## 1. Quick Decision Rules

### Should I Use a Two-Stage Design?

**YES if:**
- Covariate costs >$100-500 per subject
- You have longitudinal outcomes on all subjects before measuring covariate
- N >= 300
- Effect size >= 0.3 SD
- Can afford 15-25% sampling

**NO if:**
- Covariate is cheap or already collected
- Need >75% sampling anyway
- N < 100
- Relationship between outcome and covariate is unknown/weak

### Which Sampling Strategy?

| Condition | Recommendation |
|-----------|----------------|
| First time / regulatory | SRS (25-35%) |
| ICC >= 0.3, M <= 5 | BDS (15-25%) |
| ICC < 0.3, M >= 7 | ODS (15-25%) |
| Uncertain about effect | SRS pilot first |

### Which Analysis Method?

| Condition | Recommendation |
|-----------|----------------|
| N > 5000 | ACML (faster) |
| Prefer Bayesian | Bayesian imputation |
| Limited compute | Pathfinder or ACML |
| Default | Bayesian MCMC |

---

## 2. Recommended Parameter Settings

### ODS/BDS Design

| Parameter | Value | Notes |
|-----------|-------|-------|
| cutoff_low | 0.10 | 10th percentile |
| cutoff_high | 0.90 | 90th percentile |
| prop_low | 0.40 | 40% from low stratum |
| prop_middle | 0.20 | 20% from middle |
| prop_high | 0.40 | 40% from high stratum |
| sampling_type | "slope" | Or "intercept" if baseline effect |

### Bayesian Imputation

| Parameter | Value | Notes |
|-----------|-------|-------|
| n_chains | 4 | Minimum |
| iter_warmup | 1000-2000 | More for complex models |
| iter_sampling | 1000-2000 | More for publication |
| adapt_delta | 0.8 | Increase to 0.95 if divergences |

---

## 3. Sample Size Guidelines

For ODS/BDS with 20% sampling to detect beta_x with 80% power:

| Effect Size | Required N |
|-------------|------------|
| Small (d=0.3) | 1500-2000 |
| Medium (d=0.5) | 600-800 |
| Large (d=0.8) | 300-400 |

**Adjustments:**
- SRS instead: multiply N by 1.3-1.5
- ICC > 0.5: multiply N by 1.2
- Sampling fraction < 20%: multiply by (0.2/f)^0.7

---

## 4. Sensitivity Analyses Checklist

Before trusting results, check:

- [ ] Imputation distribution matches covariate type
- [ ] Compare cutoffs: 10th/90th vs 25th/75th
- [ ] Check MAR assumption plausibility
- [ ] MCMC diagnostics: Rhat < 1.01, ESS > 400
- [ ] Prior sensitivity (half/double prior SD)
- [ ] Asymmetry check for positive vs negative effects

---

## 5. Common Pitfalls

1. **Ignoring sampling design** - Don't use standard mixed models
2. **Too small sampling fraction** - Never below 10-15%
3. **Wrong imputation distribution** - Match to covariate type
4. **Using ODS with weak effects** - SRS better when beta_x ~ 0
5. **Ignoring ICC** - Use BDS when ICC > 0.3
6. **Not checking convergence** - Always check Rhat, ESS
7. **Overstating certainty** - Report uncertainty, not just estimates
8. **Wrong sampling variable** - Slope vs intercept matters
9. **Not planning for attrition** - Oversample by expected dropout
10. **No sensitivity analysis** - Always test robustness

---

## 6. Reporting Template

> "We used a two-stage outcome-dependent sampling design. In stage 1, we collected longitudinal outcome data on all N=1000 subjects (5 time points each). In stage 2, we measured the expensive biomarker on 200 subjects (20% sampling fraction) selected as follows: We estimated subject-specific slopes using ordinary least squares, stratified subjects into low (<10th percentile), middle (10th-90th), and high (>90th percentile) groups, and randomly sampled 40% from each tail and 20% from the middle.
>
> We analyzed the data using Bayesian imputation via the bayes2stage R package (v0.1.0), specifying a normal imputation model for the biomarker conditional on the inexpensive covariate Z. We fit the model using MCMC with 4 chains, 2000 warmup iterations, and 2000 sampling iterations. Convergence diagnostics were satisfactory (all Rhat < 1.01, minimum ESS = 450).
>
> Sensitivity analyses using (i) alternative cutoffs (25th/75th percentiles), (ii) a log-normal imputation distribution, and (iii) priors with half and double the default SD yielded qualitatively similar results (see Supplementary Table S1)."

---

## 7. Method Summary

### SRS (Simple Random Sampling)
- **Use when**: Simplicity, transparency, uncertain effects
- **Avoid when**: Budget-constrained, strong prior on effect
- **Min fraction**: 25-35%

### ODS (Outcome-Dependent Sampling)
- **Use when**: Cost efficiency critical, ICC < 0.3, M >= 7
- **Avoid when**: ICC > 0.7, weak effects, very small N
- **Min fraction**: 15-20%

### BDS (BLUP-Dependent Sampling)
- **Use when**: ICC >= 0.3, M <= 5, compute available
- **Avoid when**: M >= 10, low ICC, compute-limited
- **Min fraction**: 15-20%

### Bayesian Imputation
- **Use when**: Uncertainty quantification needed, N <= 5000
- **Avoid when**: Very large N, no compute resources

### ACML
- **Use when**: Large N, ODS design, frequentist preference
- **Avoid when**: BDS design, small N, Bayesian preference

---

## 8. When Methods Break Down

| Condition | Risk |
|-----------|------|
| N < 100 | All methods unreliable |
| sampling_fraction < 0.05 | High instability |
| \|beta_z\| > 5 | ODS bias increases |
| ICC > 0.9 | ODS loses advantage |
| gamma_1 > 0.95 | Imputation model issues |
| N_sampled < 20 | All methods fail |

---

## 9. Key Literature

- Zhou et al. (2002): ODS for continuous outcomes, Biometrics
- Rubin (1987): Multiple Imputation for Nonresponse in Surveys
- Little & Rubin (2020): Statistical Analysis with Missing Data, 3rd ed
- bayes2stage documentation: github.com/maxdrohde/bayes2stage
