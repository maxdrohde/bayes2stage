# Get parameter registry with all name mappings

Central registry mapping parameter names across different contexts: grid
names (simulation config), Stan output names, display names, and ACML
names.

## Usage

``` r
get_parameter_registry()
```

## Value

A tibble with columns: grid_name, stan_name, display_name, acml_name
