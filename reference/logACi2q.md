# Log of the Ascertainment correction piece for bivariate sampling

Calculate the log transformed ascertainment correction under a bivariate
Q_i. Also return vi

## Usage

``` r
logACi2q(yi, xi, zi, wi, beta, sigma.vc, rho.vc, sigma.e, cutpoints, SampProb)
```

## Arguments

- yi:

  n_i-response vector

- xi:

  n_i by p design matrix for fixed effects

- zi:

  n_i by q design matric for random effects (intercept and slope)

- wi:

  the pre-multiplier of yi to generate the sampling variable q_i

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

  cutpoints defining the sampling regions. (a vector of length 4 c(xlow,
  xhigh, ylow, yhigh))

- SampProb:

  Sampling probabilities from within each region (vector of length 2
  c(central region, outlying region)).

## Value

log transformed ascertainment correction
