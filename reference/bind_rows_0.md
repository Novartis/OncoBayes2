# Bind rows of multiple data frames with zero fill

A version of bind_rows out of dplyr that fills non-common columns with
zeroes instead of NA. Gives an error if any of the input data contains
NA already.

## Usage

``` r
bind_rows_0(...)
```

## Arguments

- ...:

  Data frames to combine, passed into bind_rows (see dplyr
  documentation)

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)

library(tibble)

dose_info_A <- tibble(
  group_id = "hist_A",
  drug_A = 1
)

dose_info_B <- tibble(
  group_id = "hist_B",
  drug_B = 100 * (1:2)
)

bind_rows_0(dose_info_A, dose_info_B)
#> # A tibble: 3 × 3
#>   group_id drug_A drug_B
#>   <chr>     <dbl>  <dbl>
#> 1 hist_A        1      0
#> 2 hist_B        0    100
#> 3 hist_B        0    200

## Recover user set sampling defaults
options(.user_mc_options)
```
