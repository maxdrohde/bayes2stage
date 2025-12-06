# Subject specific contribution to the lme model score (also returns marginal Vi=Cov(Y\|X))

Subject specific contribution to the lme model score (also returns
marginal Vi=Cov(Y\|X))

## Usage

``` r
li.lme.score2(subjectData, beta, sigma.vc, rho.vc, sigma.e)
```

## Arguments

- subjectData:

  a list that contains yi, xi, zi, Weights.i. Note that Weights.i is
  used for inverse probability weighting only.

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

Subject specific contribution to the log-likelihood score (also returns
marginal Vi=Cov(Y\|X))
