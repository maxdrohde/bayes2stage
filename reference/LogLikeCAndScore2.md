# Calculate the ascertainment corrected log likelihood and score

Calculate the ascertainment corrected log likelihood and score

## Usage

``` r
LogLikeCAndScore2(
  params,
  y,
  x,
  z,
  id,
  w.function,
  cutpoints,
  SampProb,
  Weights,
  ProfileCol = NA,
  Keep.liC = FALSE
)
```

## Arguments

- params:

  parameter vector c(beta, log(sigma0), log(sigma1), rho, sigmae)

- y:

  response vector

- x:

  sum(n_i) by p design matrix for fixed effects

- z:

  sum(n_i) by 2 design matric for random effects (intercept and slope)

- id:

  sum(n_i) vector of subject ids

- w.function:

  sum(n_i) vector with possible values that include "mean" (mean of
  response series), "intercept" (intercept of the regression of Yi ~ zi
  where zi is the design matrix for the random effects (solve(t.zi %*%
  zi) %*% t.zi)`[1,]`), "intercept1" (intercept of the regression of Yi
  ~ zi where zi is the design matrix for the random effects (solve(t.zi
  %*% zi) %*% t.zi)`[1,]`). "intercept2" (second intercept of the
  regression of the Yi ~ zi where zi is the design matrix for the
  bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) %*%
  t.zi)`[3,]`), "slope" (slope of the regression of Yi ~ zi where zi is
  the design matrix for the random effects (solve(t.zi %*% zi) %*%
  t.zi)`[2,]`), "slope1" (slope of the regression of Yi ~ zi where zi is
  the design matrix for the random effects (solve(t.zi %*% zi) %*%
  t.zi)`[2,]`), "slope2" (second slope of the regression of the Yi ~ zi
  where zi is the design matrix for the bivariate random effects
  (b10,b11,b20,b21) solve(t.zi %*% zi) %*% t.zi)`[4,]`) "bivar"
  (intercept and slope of the regression of Yi ~ zi where zi is the
  design matrix for the random effects (solve(t.zi %*% zi) %*%
  t.zi)`[c(1,2),]`) "mvints" (first and second intercepts of the
  bivariate regression of the Yi ~ zi where zi is the design matrix for
  the bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) %*%
  t.zi)`[c(1,3),]`) "mvslps" (first and second slopes of the bivariate
  regression of the Yi ~ zi where zi is the design matrix for the
  bivariate random effects (b10,b11,b20,b21) solve(t.zi %*% zi) %*%
  t.zi)`[c(1,3),]`). There should be one unique value per subject

- cutpoints:

  A matrix with the first dimension equal to sum(n_i). These cutpoints
  define the sampling regions (bivariate Q_i: each row is a vector of
  length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a
  vector of length 2 c(k1,k2) to define the sampling regions, i.e., low,
  middle, high). Each subject should have n_i rows of the same values.

- SampProb:

  A matrix with the first dimension equal to sum(n_i). Sampling
  probabilities from within each region (bivariate Q_i: each row is a
  vector of length 2 c(central region, outlying region); univariate Q_i:
  each row is a vector of length 3 with sampling probabilities for each
  region). Each subject should have n_i rows of the same values.

- Weights:

  Subject specific sampling weights. A vector of length sum(n_i). Not
  used unless using weighted Likelihood

- ProfileCol:

  the column number(s) for which we want fixed at the value of param.
  Maimizing the log likelihood for all other parameters while fixing
  these columns at the values of `params[ProfileCol]`

- Keep.liC:

  If TRUE outputs subject specific conditional log lileihoods to be used
  for the imputation procedure described in the AOAS paper keep z
  sum(n_i) by 2 design matric for random effects (intercept and slope)

## Value

The conditional log likelihood with a "gradient" attribute (if
Keep.liC=FALSE) and subject specific contributions to the conditional
likelihood if Keep.liC=TRUE).
