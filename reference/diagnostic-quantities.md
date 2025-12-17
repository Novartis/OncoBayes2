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
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)

example_model("single_agent", silent = TRUE)
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

head(log_posterior(blrmfit))
#>   Chain Iteration     Value
#> 1     1         1 -17.00453
#> 2     1         2 -15.89877
#> 3     1         3 -13.37728
#> 4     1         4 -13.70039
#> 5     1         5 -15.04728
#> 6     1         6 -18.13408

np <- nuts_params(blrmfit)
str(np)
#> 'data.frame':    60 obs. of  4 variables:
#>  $ Chain    : int  1 1 1 1 1 1 1 1 1 1 ...
#>  $ Iteration: int  1 2 3 4 5 6 7 8 9 10 ...
#>  $ Parameter: Factor w/ 6 levels "accept_stat__",..: 1 1 1 1 1 1 1 1 1 1 ...
#>  $ Value    : num  0.999 1 0.999 0.997 1 ...
# extract the number of divergence transitions
sum(subset(np, Parameter == "divergent__")$Value)
#> [1] 0

head(rhat(blrmfit))
#>        mu_log_beta[log(drug_A/dref),intercept] 
#>                                       1.391186 
#>        mu_log_beta[log(drug_A/dref),log_slope] 
#>                                       1.384032 
#>     tau_log_beta[1,log(drug_A/dref),intercept] 
#>                                            NaN 
#>     tau_log_beta[1,log(drug_A/dref),log_slope] 
#>                                            NaN 
#> beta_group[trial_A,log(drug_A/dref),intercept] 
#>                                       1.391186 
#>     beta_group[trial_A,log(drug_A/dref),slope] 
#>                                       1.238610 
head(neff_ratio(blrmfit))
#>        mu_log_beta[log(drug_A/dref),intercept] 
#>                                      0.6793941 
#>        mu_log_beta[log(drug_A/dref),log_slope] 
#>                                      0.5941920 
#>     tau_log_beta[1,log(drug_A/dref),intercept] 
#>                                            NaN 
#>     tau_log_beta[1,log(drug_A/dref),log_slope] 
#>                                            NaN 
#> beta_group[trial_A,log(drug_A/dref),intercept] 
#>                                      0.6793941 
#>     beta_group[trial_A,log(drug_A/dref),slope] 
#>                                      0.8094134 

## Recover user set sampling defaults
options(.user_mc_options)
```
