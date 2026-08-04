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
.user_mc_options <- options()

example_model("single_agent", silent = TRUE)

posterior_interval(blrmfit)
#>                                                      2.5%     97.5%
#> mu_log_beta[log(drug_A/dref),intercept]        -1.4601463 4.3873983
#> mu_log_beta[log(drug_A/dref),log_slope]        -0.6559907 1.6527162
#> tau_log_beta[1,log(drug_A/dref),intercept]      0.0000000 0.0000000
#> tau_log_beta[1,log(drug_A/dref),log_slope]      0.0000000 0.0000000
#> rho_log_beta[log(drug_A/dref)]                 -0.9474960 0.9564961
#> beta_group[trial_A,log(drug_A/dref),intercept] -1.4601463 4.3873983
#> beta_group[trial_A,log(drug_A/dref),slope]      0.5189278 5.2211429

## Recover user set sampling defaults
options(.user_mc_options)
```
