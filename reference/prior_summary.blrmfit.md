# Summarise model prior

Extracts a summary of the prior in a structured data format.

## Usage

``` r
# S3 method for class 'blrmfit'
prior_summary(object, digits = 2, ...)
```

## Arguments

- object:

  `blrmfit` (`blrm_trial`) object as returned from
  [`blrm_exnex()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_exnex.md)
  ([`blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_trial.md))
  analysis

- digits:

  number of digits to show

- ...:

  ignored by the function

## Value

Returns an analysis specific list, which has it's own `print` function.
The returned list contains arrays which represent the prior in a
structured format.

## Details

The summary of the prior creates a structured representation of the
specified prior from a
[`blrm_exnex()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_exnex.md)
([`blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_trial.md))
analysis.

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)

## run combo2 analysis which defines blrmfit model object
example_model("combo2", silent = TRUE)
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess

prior_summary(blrmfit)
#> Bayesian Logistic Regression Model with EXchangeability-NonEXchangeability
#> 
#> Mixture configuration
#> ---------------------
#> EXNEX components : 0 
#> component
#> I(log(drug_A/dref[1])) I(log(drug_B/dref[2])) 
#>                      0                      0 
#> 
#> EXNEX interactions: 0 
#> interaction
#> I(drug_A/dref[1] * drug_B/dref[2]) 
#>                                  0 
#> 
#> Prior probability for exchangeability per group
#>           component
#> group      I(log(drug_A/dref[1])) I(log(drug_B/dref[2]))
#>   trial_A                       1                      1
#>   trial_B                       1                      1
#>   IIT                           1                      1
#>   trial_AB                      1                      1
#> 
#>           interaction
#> group      I(drug_A/dref[1] * drug_B/dref[2])
#>   trial_A                                   1
#>   trial_B                                   1
#>   IIT                                       1
#>   trial_AB                                  1
#> 
#> EXchangable hyperparameter priors
#> ---------------------------------
#> Component parameters
#> Mean mu_log_beta
#>                               prior weight m_intercept m_log_slope s_intercept s_log_slope  rho
#> component              mix                                                                     
#> I(log(drug_A/dref[1])) comp_1          1.0        -1.4         0.0         2.0         1.0  0.0
#> I(log(drug_B/dref[2])) comp_1          1.0        -1.4         0.0         2.0         1.0  0.0
#> 
#> Heterogeneity tau_log_beta (log-normal)
#>                                         prior weight m_tau_intercept m_tau_log_slope s_tau_intercept s_tau_log_slope   rho
#> stratum   component              mix                                                                                      
#> stratum_1 I(log(drug_A/dref[1])) comp_1         1.00           -1.39           -2.08            0.71            0.71  0.00
#>           I(log(drug_B/dref[2])) comp_1         1.00           -1.39           -2.08            0.71            0.71  0.00
#> 
#> Correlation LKJ
#> component
#> I(log(drug_A/dref[1])) I(log(drug_B/dref[2])) 
#>                      1                      1 
#> 
#> Interaction parameters
#> Mean mu_eta
#>       prior   w m[1] s[1]
#> mix                      
#> comp1       1.0  0.0  1.1
#> 
#> Heterogeneity tau_eta (log-normal)
#>                 prior     w  m[1]  s[1]
#> stratum   mix                          
#> stratum_1 comp1        1.00 -2.08  0.71
#> 
#> Correlation LKJ
#> interaction 
#>           1 
#> 
#> NonEXchangable priors
#> ---------------------
#> Component parameters
#> Mean mu_log_beta
#>                               prior weight m_intercept m_log_slope s_intercept s_log_slope  rho
#> component              mix                                                                     
#> I(log(drug_A/dref[1])) comp_1          1.0        -1.4         0.0         2.0         1.0  0.0
#> I(log(drug_B/dref[2])) comp_1          1.0        -1.4         0.0         2.0         1.0  0.0
#> 
#> Interaction parameters
#> Mean mu_eta
#>       prior   w m[1] s[1]
#> mix                      
#> comp1       1.0  0.0  1.1

prior_sum <- prior_summary(blrmfit)
names(prior_sum)
#>  [1] "has_inter"         "is_EXNEX_comp"     "EX_prob_comp"     
#>  [4] "EX_mu_log_beta"    "NEX_mu_log_beta"   "EX_tau_log_beta"  
#>  [7] "EX_corr_eta_comp"  "tau_dist"          "is_EXNEX_inter"   
#> [10] "EX_prob_inter"     "EX_mu_eta"         "NEX_mu_eta"       
#> [13] "EX_tau_eta"        "EX_corr_eta_inter" "num_strata"       
#> [16] "num_groups"       

## the entries of the prior list are labelled arrays
dimnames(prior_sum$EX_mu_log_beta)
#> $prior
#> [1] "weight"      "m_intercept" "m_log_slope" "s_intercept" "s_log_slope"
#> [6] "rho"        
#> 
#> $component
#> [1] "I(log(drug_A/dref[1]))" "I(log(drug_B/dref[2]))"
#> 
#> $mix
#> [1] "comp_1"
#> 

## Recover user set sampling defaults
options(.user_mc_options)
```
