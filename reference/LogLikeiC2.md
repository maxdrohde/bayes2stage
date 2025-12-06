# Calculate the ss contributions to the conditional likelihood for the univariate and bivariate sampling cases.

Calculate the ss contributions to the conditional likelihood for the
univariate and bivariate sampling cases.

## Usage

``` r
LogLikeiC2(subjectData, beta, sigma.vc, rho.vc, sigma.e)
```

## Arguments

- subjectData:

  a list containing: yi, xi, zi, Weights.i, w.function.i, SampProb.i,
  cutpoints.i

- beta:

  mean model parameter p-vector

- sigma.vc:

  vector of variance components on standard deviation scale

- rho.vc:

  vector of correlations among the random effects. The length should be
  q choose 2

- sigma.e:

  std dev of the measurement error distribution

## Value

ss contributions to the conditional log likelihood. This is an internal
function used by LogLikeC2
