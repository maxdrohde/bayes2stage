# Ascertainment correction piece for univariate sampling

Calculate the (not yet log transformed) ascertainment correction under a
univariate Q_i

## Usage

``` r
ACi1q(cutpoints, SampProb, mu_q, sigma_q)
```

## Arguments

- cutpoints:

  cutpoints defining the sampling regions. (a vector of length 2)

- SampProb:

  Sampling probabilities from within each region (vector of length 3).

- mu_q:

  a scalar for the mean value of the Q_i distribution

- sigma_q:

  a scalar for the standard deviation of the Q_i distribution

## Value

Not yet log transformed ascertainment correction
