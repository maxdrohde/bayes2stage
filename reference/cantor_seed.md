# Generate a unique seed from two indices using Cantor pairing

Uses the Cantor pairing function to generate a unique integer seed from
two positive integer indices. Useful for reproducible simulations where
each (i, j) combination needs a unique seed.

## Usage

``` r
cantor_seed(i, j)
```

## Arguments

- i:

  First positive integer index

- j:

  Second positive integer index

## Value

A unique integer seed

## Examples

``` r
cantor_seed(1, 1)
#> [1] 4
cantor_seed(5, 10)
#> [1] 130
```
