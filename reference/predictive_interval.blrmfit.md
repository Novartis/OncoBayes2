# Posterior predictive intervals

Posterior predictive intervals of the model.

## Usage

``` r
# S3 method for class 'blrmfit'
predictive_interval(object, prob = 0.95, newdata, ...)
```

## Arguments

- object:

  fitted model object

- prob:

  central probability mass to report, i.e. the quantiles 0.5-prob/2 and
  0.5+prob/2 are displayed. Multiple central widths can be specified.

- newdata:

  optional data frame specifying for what to predict; if missing, then
  the data of the input model `object` is used

- ...:

  not used in this function

## Value

Matrix with as many rows as the input data set and two columns which
contain the lower and upper quantile corresponding to the central
probability mass `prob` for the number of responses of the predictive
distribution.

## Details

Reports for each row of the input data set the predictive interval
according to the fitted model.

## Examples

``` r
.user_mc_options <- options()

example_model("single_agent", silent = TRUE)

predictive_interval(blrmfit)
#>   2.5% 97.5%
#> 1    0     0
#> 2    0     1
#> 3    0     2
#> 4    0     2
#> 5    0     2

## Recover user set sampling defaults
options(.user_mc_options)
```
