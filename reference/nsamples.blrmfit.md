# Return the number of posterior samples

Return the number of posterior samples

## Usage

``` r
# S3 method for class 'blrmfit'
nsamples(object, ...)
```

## Arguments

- object:

  fitted model object

- ...:

  not used in this function

## Examples

``` r
.user_mc_options <- options()


## run single-agent analysis which defines blrmfit model object
example_model("single_agent", silent = TRUE)

nsamples(blrmfit)
#> [1] 4000

## Recover user set sampling defaults
options(.user_mc_options)
```
