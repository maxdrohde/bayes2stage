# Calculate V_i = Z_i D t(Z_i) + sig_e^2 `I_(n_i)`

Calculate V_i = Z_i D t(Z_i) + sig_e^2 `I_(n_i)`

## Usage

``` r
vi.calc(zi, sigma.vc, rho.vc, sigma.e)
```

## Arguments

- zi:

  n_i by q design matrix for the random effects

- sigma.vc:

  vector of variance components on standard deviation scale

- rho.vc:

  vector of correlations among the random effects. The length should be
  q choose 2

- sigma.e:

  std dev of the measurement error distribution

## Value

V_i
