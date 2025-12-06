# Calculate a subject-specific contribution to a log-likelihood for longitudinal normal data

Calculate a subject-specific contribution to a log-likelihood for
longitudinal normal data

## Usage

``` r
li.lme(yi, xi, beta, vi)
```

## Arguments

- yi:

  n_i-response vector

- xi:

  n_i by p design matrix for fixed effects

- beta:

  mean model parameter vector

- vi:

  the variance covariance matrix (ZDZ+Sige2\*I)

## Value

subject specific contribution to the log-likelihood
