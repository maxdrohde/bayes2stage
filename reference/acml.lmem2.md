# Fitting function: ACML or WL for a linear mixed effects model (random intercept and slope)

Fitting function: ACML or WL for a linear mixed effects model (random
intercept and slope)

## Usage

``` r
acml.lmem2(
  formula.fixed,
  formula.random,
  data,
  id,
  w.function = "mean",
  InitVals,
  cutpoints,
  SampProb,
  Weights,
  ProfileCol = NA
)
```

## Arguments

- formula.fixed:

  formula for the fixed effects (of the form y~x)

- formula.random:

  formula for the random effects (of the form ~z). Right now this model
  only fits random intercept and slope models.

- data:

  data frame that should contain everything in formula.fixed,
  formula.random, id, and Weights. It does not include: w.function,
  cutpoints, SampProb

- id:

  sum(n_i) vector of subject ids (a variable contained in data)

- w.function:

  sum(n_i) vector with possible values that include "mean" (mean of
  response series), "intercept" (intercept of the regression of Yi ~ zi
  where zi is the design matrix for the random effects (solve(t.zi\* zi)
  \* t.zi)`[1,]`), "intercept1" (intercept of the regression of Yi ~ zi
  where zi is the design matrix for the random effects (solve(t.zi \*
  zi) \* t.zi)`[1,]`). "intercept2" (second intercept of the regression
  of the Yi ~

- InitVals:

  starting values for c(beta, log(sigma0), log(sigma1),
  log((1+rho)/(1-rho)), log(sigmae))

- cutpoints:

  A matrix with the first dimension equal to sum(n_i). These cutpoints
  define the sampling regions for individual subjects. If using a low,
  medium, high, sampling scheme, this is a sum(n_i) by 2 matrix that
  must be a distinct object not contained in the dat dataframe. Each row
  is a vector of length 2 c(k1,k2) to define the sampling regions, i.e.,
  low, middle, high. If using a square doughnut design this should be
  sum(n_i) by 4 matrix (var1lower, var1upper, var2lower, var2upper).
  Each subject should have n_i rows of the same values.

- SampProb:

  A matrix with the first dimension equal to sum(n_i). Sampling
  probabilities from within each region. For low medium high sampling,
  each row is a vector of length 3 with sampling probabilities for each
  region. For bivariate stratum sampling each row is a vector of length
  2 with sampling probabilities for the inner and outer strata. Each
  subject should have n_i rows of the same values. Not in data.

- Weights:

  Subject specific sampling weights. A vector of length sum(n_i). This
  should be a variable in the data dataframe. It should only be used if
  doing IPWL. Note if doing IPWL, only use robcov (robust variances) and
  not covar. If not doing IPWL, this must be a vectors of 1s.

- ProfileCol:

  the column number(s) for which we want fixed at the value of param.
  Maximizing the log likelihood for all other parameters while fixing
  these columns at the values of `InitVals[ProfileCol]`

## Value

Ascertainment corrected Maximum likelihood: Ests, covar, LogL, code,
robcov
