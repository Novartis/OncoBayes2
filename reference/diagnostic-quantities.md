# Extract Diagnostic Quantities of OncoBayes2 Models

Extract quantities that can be used to diagnose sampling behavior of the
algorithms applied by Stan at the back-end of OncoBayes2.

## Usage

``` r
# S3 method for class 'blrmfit'
log_posterior(object, ...)

# S3 method for class 'blrmfit'
nuts_params(object, pars = NULL, ...)

# S3 method for class 'blrmfit'
rhat(object, pars = NULL, ...)

# S3 method for class 'blrmfit'
neff_ratio(object, pars = NULL, ...)
```

## Arguments

- object:

  A `blrmfit` or `blrmtrial` object.

- ...:

  Arguments passed to individual methods.

- pars:

  An optional character vector of parameter names. For `nuts_params`
  these will be NUTS sampler parameter names rather than model
  parameters. If `pars` is omitted all parameters are included.

## Value

The exact form of the output depends on the method.

## Details

For more details see
[`bayesplot::bayesplot-extractors()`](https://mc-stan.org/bayesplot/reference/bayesplot-extractors.html).

## Examples

``` r
.user_mc_options <- options()

example_model("single_agent", silent = TRUE)

head(log_posterior(blrmfit))
#>   Chain Iteration     Value
#> 1     1         1 -13.30890
#> 2     1         2 -14.46850
#> 3     1         3 -14.05445
#> 4     1         4 -13.82799
#> 5     1         5 -15.67133
#> 6     1         6 -15.40357

np <- nuts_params(blrmfit)
str(np)
#> 'data.frame':    24000 obs. of  4 variables:
#>  $ Chain    : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ Iteration: int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ Parameter: Factor w/ 6 levels "accept_stat__",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ Value    : num  1 0.988 0.999 0.997 0.989 ...
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value)
#> [1] 0

head(rhat(blrmfit))
#>                     log_beta_raw[1,1,1]                     log_beta_raw[2,1,1] 
#>                               1.0002393                               1.0022189 
#>                     log_beta_raw[1,1,2]                     log_beta_raw[2,1,2] 
#>                               0.9996421                               1.0021369 
#> mu_log_beta[log(drug_A/dref),intercept] mu_log_beta[log(drug_A/dref),log_slope] 
#>                               1.0013644                               1.0011005 
head(neff_ratio(blrmfit))
#>                     log_beta_raw[1,1,1]                     log_beta_raw[2,1,1] 
#>                               1.2565507                               1.1143313 
#>                     log_beta_raw[1,1,2]                     log_beta_raw[2,1,2] 
#>                               1.1332693                               1.0414032 
#> mu_log_beta[log(drug_A/dref),intercept] mu_log_beta[log(drug_A/dref),log_slope] 
#>                               0.6854192                               0.7189484 

## Recover user set sampling defaults
options(.user_mc_options)
```
