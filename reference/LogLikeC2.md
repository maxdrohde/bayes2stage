# Calculate the conditional likelihood for the univariate and bivariate sampling cases across all subjects (Keep.liC=FALSE) or the subject specific contributions to the conditional likelihood along with the log-transformed ascertainment correction for multiple imputation (Keep.liC=TRUE).

Calculate the conditional likelihood for the univariate and bivariate
sampling cases across all subjects (Keep.liC=FALSE) or the subject
specific contributions to the conditional likelihood along with the
log-transformed ascertainment correction for multiple imputation
(Keep.liC=TRUE).

## Usage

``` r
LogLikeC2(
  y,
  x,
  z,
  w.function,
  id,
  beta,
  sigma.vc,
  rho.vc,
  sigma.e,
  cutpoints,
  SampProb,
  Weights,
  Keep.liC = FALSE
)
```

## Arguments

- y:

  response vector

- x:

  sum(n_i) by p design matrix for fixed effects

- z:

  sum(n_i) by q design matrix for random effects

- w.function:

  sum(n_i) vector with possible values that include "mean" "intercept"
  "slope" and "bivar." There should be one unique value per subject

- id:

  sum(n_i) vector of subject ids

- beta:

  mean model parameter p-vector

- sigma.vc:

  vector of variance components on standard deviation scale

- rho.vc:

  vector of correlations among the random effects. The length should be
  q choose 2

- sigma.e:

  std dev of the measurement error distribution

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

- Keep.liC:

  If FALSE, the function returns the conditional log likelihood across
  all subjects. If TRUE, subject specific contributions and
  exponentiated subject specific ascertainment corrections are returned
  in a list.

## Value

If Keep.liC=FALSE, conditional log likelihood. If Keep.liC=TRUE, a
two-element list that contains subject specific likelihood contributions
and exponentiated ascertainment corrections.
