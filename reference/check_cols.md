# Check for required columns in a data frame

Validates that a data frame contains all required columns.

## Usage

``` r
check_cols(data, required_cols)
```

## Arguments

- data:

  A data frame to check

- required_cols:

  A character vector of required column names

## Value

Invisibly returns TRUE if all columns are present; otherwise throws an
error
