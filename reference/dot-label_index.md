# Utility function to label parameter indices according to factor levels.

Utility function to label parameter indices according to factor levels.

## Usage

``` r
.label_index(stanfit, par, ...)
```

## Arguments

- stanfit:

  stan fit which names are being modified

- par:

  parameter selected

- ...:

  must include as many factors as there are indices which are used in
  the order given to translate indices to text labels
