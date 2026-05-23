# Translate parameter names between contexts

Translate parameter names between contexts

## Usage

``` r
translate_parameter_names(names, from = "stan", to = "display")
```

## Arguments

- names:

  Character vector of parameter names to translate

- from:

  Source context: "grid", "stan", "display", or "acml"

- to:

  Target context: "grid", "stan", "display", or "acml"

## Value

Character vector of translated names (preserves original if no mapping
found)
