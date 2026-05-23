# Extract summary statistics from model fit objects

S3 generic for extracting summary statistics from various fit object
types. Supports CmdStanR fit types and ACML fits.

## Usage

``` r
model_summary(object, ...)

# S3 method for class 'CmdStanMCMC'
model_summary(object, ...)

# S3 method for class 'CmdStanPathfinder'
model_summary(object, ...)

# S3 method for class 'CmdStanLaplace'
model_summary(object, ...)

# S3 method for class 'CmdStanMLE'
model_summary(object, ...)

# S3 method for class 'acml_fit'
model_summary(object, ...)

# Default S3 method
model_summary(object, ...)
```

## Arguments

- object:

  A model fit object

- ...:

  Additional arguments passed to methods

## Value

A data frame with summary statistics for each parameter
