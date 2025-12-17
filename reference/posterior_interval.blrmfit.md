# Posterior intervals

Posterior intervals of all model parameters.

## Usage

``` r
# S3 method for class 'blrmfit'
posterior_interval(object, prob = 0.95, ...)
```

## Arguments

- object:

  fitted model object

- prob:

  central probability mass to report, i.e. the quantiles 0.5-prob/2 and
  0.5+prob/2 are displayed. Multiple central widths can be specified.

- ...:

  not used in this function

## Value

Matrix of two columns for the central probability interval `prob` for
all parameters of the model.

## Details

Reports the quantiles of posterior parameters which correspond to the
central probability mass specified. The output includes the posterior of
the hyper-parameters and the posterior of each group estimate.

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

posterior_interval(blrmfit)
#>                                                      2.5%     97.5%
#> mu_log_beta[log(drug_A/dref),intercept]        -0.8640091 3.8839930
#> mu_log_beta[log(drug_A/dref),log_slope]        -0.4132765 1.5579169
#> tau_log_beta[1,log(drug_A/dref),intercept]      0.0000000 0.0000000
#> tau_log_beta[1,log(drug_A/dref),log_slope]      0.0000000 0.0000000
#> rho_log_beta[log(drug_A/dref)]                 -0.8234919 0.8211436
#> beta_group[trial_A,log(drug_A/dref),intercept] -0.8640091 3.8839930
#> beta_group[trial_A,log(drug_A/dref),slope]      0.7324769 4.7797538

## Recover user set sampling defaults
options(.user_mc_options)
```
