# Statistical Theory Behind Simulation Hypotheses

## H1: Efficiency Gains Scale with Effect Size

### Mathematical Basis

The efficiency of outcome-dependent sampling (ODS) relative to balanced or simple random sampling stems from **Fisher information theory**. For a parameter beta in a regression model, the Fisher Information is:

I(beta) = E[-d^2 l / d beta^2]

where l is the log-likelihood. The asymptotic variance of the MLE is Var(beta_hat) ~ 1/I(beta).

In ODS designs, sampling is concentrated in regions where the **score function** dl/d beta has high variance. For a linear mixed model:

y_ij = beta_x * x_i + beta_z * z_i + beta_t * t_ij + b_i + epsilon_ij

the contribution to Fisher information from subject i is approximately proportional to:

I_i(beta_z) ~ Var(dl_i/d beta_z | data) ~ z_i^2 * Var(y_ij | z_i, x_i)

**When |beta_z| is large**: Subjects with extreme z_i have systematically extreme outcomes y_ij. ODS, by sampling extreme y values, automatically enriches the sample with subjects having extreme z_i values. These subjects contribute disproportionately to I(beta_z).

**When beta_z ~ 0**: Extreme y values occur for reasons unrelated to z. ODS provides no informational advantage.

### Theoretical Framework

This relates to **Neyman optimal allocation** in survey sampling theory (Neyman 1934). For estimating a population mean, sample sizes should be proportional to stratum variances:

n_h ~ N_h * sigma_h

Zhou et al. (2002) showed that sampling from the tails of the outcome distribution is asymptotically optimal for estimating exposure effects.

### Precise Statistical Quantities

**Primary comparison**:
- Asymptotic relative efficiency (ARE) = Var(beta_z_hat^BDS) / Var(beta_z_hat^ODS)
- Empirical relative efficiency = (SE_BDS / SE_ODS)^2

**Expected relationship**:
ARE(beta_z) = 1 + c * |beta_z|^alpha

where c > 0 and alpha ~ 1 or 2.

### Required Assumptions

1. Correct model specification
2. Conditional independence given covariates
3. Sufficient overlap (positivity)
4. No unmeasured confounding
5. Large sample for asymptotic properties

---

## H2: Asymmetry Around Zero

### Mathematical Basis

In the model with beta_x = -1, beta_t = -1, and beta_z varying, the **conditional distribution** given z has asymmetric tails when beta_z != 0.

**Lower tail (small y)**: When beta_z = -4, subjects with large z_i contribute to lower tail. When beta_z = +4, subjects with small z_i contribute.

Because x and z may be correlated through gamma_1 in the imputation model, the **joint distribution** of covariates differs between tails. This creates **selection-induced correlation**.

### Theoretical Framework

This relates to **truncated bivariate normal theory**. When sampling is based on y > c (upper tail):

f(x, z | y > c) ~ f(x, z, y) * I(y > c)

which induces correlation between x and z even if marginally independent (Heckman 1979).

### Expected Relationship

If gamma_1 != 0 and beta_x != 0:
Bias(beta_z_hat | beta_z = -k) != Bias(beta_z_hat | beta_z = +k)

The magnitude of asymmetry scales with gamma_1 * beta_x.

---

## H3: Null Effect Behavior

### Mathematical Basis

When beta_z = 0, z is **causally irrelevant** for y. ODS samples based on S = I(y in extreme tails), where:

P(S = 1 | z) = P(y in tails | z) = P(y in tails) when beta_z = 0

This means sampling is **independent of z**, satisfying **missing at random (MAR)**.

### Theoretical Framework

Under MAR, likelihood-based methods are **unbiased** (Rubin 1976). The likelihood factors as:

L(theta | observed) ~ L_data(theta) * L_sampling(phi)

The sampling mechanism is **ignorable** for likelihood inference.

**However**, precision differs. ODS samples extreme y values due to extreme x, t, b_i, or epsilon. This may **reduce** Var(z | sampled) compared to SRS, increasing variance of beta_z estimates.

---

## H4: Bias-Variance Tradeoff Across Effect Sizes

### Mathematical Basis

ODS sampling creates a **collider stratification** problem:

```
x -> y <- z
     ^
     S (sampling indicator)
```

Conditioning on S induces a **spurious association** between x and z (Pearl 2009, Hernan et al. 2004).

The induced bias can be approximated as:
Bias(beta_z_hat) ~ beta_x * Cov(x, z | y in tails, observed covariates)

When beta_z is large, collider bias is **diluted** because the signal (z -> y) dominates.

**However**, if the imputation model correctly specifies the x-z relationship, Bayesian imputation can **remove** this bias.

### Bias Bound (Ding & Miratrix 2015)

|Bias(beta_z_hat)| <= C * sqrt(1/n_lower + 1/n_upper)

---

## H5: Imputation Model Strength

### Mathematical Basis

In the imputation model x_i ~ N(gamma_0 + gamma_1 * z_i, sigma_x^2), the R^2 is:

R^2 = gamma_1^2 * Var(z) / [gamma_1^2 * Var(z) + sigma_x^2]

**When gamma_1 = 0**: z provides no information about x, so imputation adds **no information**.

**When gamma_1 is large**: z strongly predicts x, so imputation **recovers information**.

The **effective sample size** for estimating beta_x:

n_eff = n_observed + n_imputed * R^2

### Multiple Imputation Efficiency (Rubin 1987)

Var(beta_hat | imputation) = Var(beta_hat | full data) * (1 + lambda/m)

where lambda = (Var_between / Var_within) is the fraction of missing information.

---

## H6: Sampling Fraction Effects

### Mathematical Basis

With sampling fraction f, expected information under SRS is:
E[I_SRS] = f * I_total

Under ODS:
E[I_ODS] = f * E[I_i | y_i in tails]

**When f is small** (e.g., 0.1): ODS selects the 10% with highest I_i, so E[I_ODS] >> E[I_SRS].

**When f is large** (e.g., 0.75): Difference only in the 25% not sampled, so E[I_ODS] - E[I_SRS] -> 0 as f -> 1.

### Design Effect

Deff(f) = Var_ODS(f) / Var_SRS(f)

- Deff(f) << 1 when f is small (ODS much better)
- Deff(f) -> 1 as f -> 1 (designs converge)

---

## H7: Random Effects Magnitude

### Mathematical Basis

With random intercepts b_i ~ N(0, sigma_b^2), the ICC is:

ICC = sigma_b^2 / (sigma_b^2 + sigma_epsilon^2)

**When ICC is high**: Between-subject variation dominates. Extreme y values primarily reflect extreme b_i, not extreme covariate values.

ODS essentially samples subjects with extreme b_i, providing:
- High information about sigma_b
- Low information about beta_z

**Effective sample size in clustered data** (Eldridge et al. 2006):

n_eff = n_total / [1 + (m - 1) * ICC]

**Design effect for ODS under high ICC**:
Deff ~ 1 + c * ICC

where c > 0, showing ODS loses advantage as ICC increases.

---

## H8: Sample Size Scaling

### Mathematical Basis

By asymptotic theory:
beta_hat ~ N(beta, I(beta)^{-1}/n)

For ODS and SRS:
SE_ODS ~ 1/sqrt(n * I_ODS)
SE_SRS ~ 1/sqrt(n * I_SRS)

The **ratio** SE_ODS / SE_SRS = sqrt(I_SRS / I_ODS) is **independent of n** asymptotically.

At small n, finite-sample bias may differ, but as n -> infinity:
- Bias -> 0 for both (consistency)
- SE ratio -> constant

---

## H9: ODS Cutoff Sensitivity

### Mathematical Basis

ODS samples from y < q_low and y > q_high.

**More extreme cutoffs** (e.g., 2.5th/97.5th):
- Concentrate on most extreme z values (if beta_z large)
- Increase Fisher information per subject
- Increase **truncation severity**, amplifying collider bias

The bias in truncated samples (Heckman 1979):
Bias(beta_hat) ~ sigma_ez * lambda(cutoff)

where lambda(c) = phi(c) / Phi(c) is the **inverse Mills ratio**.

**As cutoffs become more extreme**:
- lambda increases -> Bias increases
- Variance decreases (more informative observations)

This creates a **bias-variance tradeoff**.

---

## H10: X Distribution Effects

### Mathematical Basis

**Continuous x** (normal): Imputation uses linear regression with Gaussian posterior.

**Binary x** (Bernoulli): Imputation uses logistic regression with Bernoulli posterior.

Key difference: Binary imputation introduces non-linearity and discreteness.

**Mutual information**:
I(z; x) = H(x) - H(x | z)

For continuous x: H(x) = 0.5 * log(2*pi*e * sigma_x^2)
For binary x: H(x) = -p*log(p) - (1-p)*log(1-p) <= 1 bit

Binary x provides **less information** to recover from z, so efficiency gains from imputation are reduced.

---

## Key Literature

- **Zhou et al. (2002)**: ODS for continuous outcomes (Biometrics)
- **Rubin (1976, 1987)**: Missing data theory
- **Pearl (2009)**: Causality and collider bias
- **Hernan et al. (2004)**: Selection bias (Epidemiology)
- **Heckman (1979)**: Selection models (Econometrica)
- **Eldridge et al. (2006)**: Effective sample size in clusters
- **Ding & Miratrix (2015)**: Sensitivity analysis for M-bias
