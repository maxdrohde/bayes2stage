# Create a list of subject-specific data

Create a list of subject-specific data

## Usage

``` r
CreateSubjectData(id, y, x, z, Weights, SampProb, cutpoints, w.function)
```

## Arguments

- id:

  sum(n_i) vector of subject ids

- y:

  response vector

- x:

  sum(n_i) by p design matrix for fixed effects

- z:

  sum(n_i) by 2 design matric for random effects (intercept and slope)

- Weights:

  Subject specific sampling weights. A vector of length sum(n_i). Not
  used unless using weighted Likelihood

- SampProb:

  A matrix with the first dimension equal to sum(n_i). Sampling
  probabilities from within each region (bivariate Q_i: each row is a
  vector of length 2 c(central region, outlying region); univariate Q_i:
  each row is a vector of length 3 with sampling probabilities for each
  region). Each subject should have n_i rows of the same values.

- cutpoints:

  A matrix with the first dimension equal to sum(n_i). These cutpoints
  define the sampling regions (bivariate Q_i: each row is a vector of
  length 4 c(xlow, xhigh, ylow, yhigh); univariate Q_i: each row is a
  vector of length 2 c(k1,k2) to define the sampling regions, i.e., low,
  middle, high). Each subject should have n_i rows of the same values.

- w.function:

  sum(n_i) vector with possible values that include "mean" "intercept"
  "slope" and "bivar." There should be one unique value per subject
