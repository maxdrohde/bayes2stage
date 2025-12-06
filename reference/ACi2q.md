# Ascertainment correction piece for bivariate sampling

Calculate the (not yet log transformed) ascertainment correction under a
bivariate Q_i

## Usage

``` r
ACi2q(cutpoints, SampProb, mu_q, sigma_q)
```

## Arguments

- cutpoints:

  cutpoints defining the sampling regions. (a vector of length 4:
  c(xlow, xhigh, ylow, yhigh))

- SampProb:

  Sampling probabilities from within each of two sampling regions;
  central region and outlying region (vector of length 2).

- mu_q:

  a 2-vector for the mean value of the bivariate Q_i distribution.

- sigma_q:

  a 2 by 2 covariance matrix for the bivariate Q_i distribution.

## Value

Not yet log transformed ascertainment correction
