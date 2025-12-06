# Gradient of the log of the ascertainment correction piece for sampling based on univariate Q_i

Calculate the gradient of the log transformed ascertainment correction
for sampling based on univariate Q_i

## Usage

``` r
logACi1q.score2(subjectData, beta, sigma.vc, rho.vc, sigma.e)
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

gradient of the log transformed ascertainment correction under
univariate \$Q_i\$
