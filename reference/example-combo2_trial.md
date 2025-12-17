# Two-drug combination example using BLRM Trial

Example using
[`blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_trial.md)
to guide the built-in two-drug combination study example.

## Details

[`blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_trial.md)
is used to collect and store all relevant design information for the
example. Subsequent use of the
[`update.blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/update.blrm_trial.md)
command allows convenient model fitting via
[`blrm_exnex()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_exnex.md).
The
[`summary.blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/summary.blrm_trial.md)
method allows exploration of the design and modeling results.

To run this example, use `example_model("combo2_trial")`. See
[`example_model()`](https://opensource.nibr.com/OncoBayes2/reference/example_model.md).

## See also

Other blrm_trial combo2 example:
[`blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_trial.md),
[`dose_info_combo2`](https://opensource.nibr.com/OncoBayes2/reference/dose_info_combo2.md),
[`drug_info_combo2`](https://opensource.nibr.com/OncoBayes2/reference/drug_info_combo2.md)

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)

library(tibble)
library(dplyr)
library(tidyr)

# Combo2 example using blrm_trial functionality

# construct initial blrm_trial object from built-in example datasets
combo2_trial_setup <- blrm_trial(
  data = hist_combo2,
  dose_info = dose_info_combo2,
  drug_info = drug_info_combo2,
  simplified_prior = FALSE
)
#> No stratum defined - assigning all groups to single stratum "all"
#> Please configure blrm_exnex using the update() function.

# summary of dimensionality of data structures
dims <- summary(combo2_trial_setup, "dimensionality")

# Fit the initial model with the historical data and fully specified prior


combo2_trial_start <- update(
   combo2_trial_setup,
   ## bivariate normal prior for drug A and drug B of intercept and
   ## log-slope
   prior_EX_mu_comp =
     replicate(2,
               mixmvnorm(c(1,
                           logit(0.2), 0,
                           diag(c(2^2, 1))))
             , FALSE),
   prior_EX_tau_comp =
     replicate(2,
               mixmvnorm(c(1,
                           log(0.25), log(0.125),
                           diag(c(log(4)/1.96, log(4)/1.96)^2)))
             , FALSE),
   prior_EX_mu_inter = mixmvnorm(c(1, 0, 1.121^2)),
   prior_EX_tau_inter = mixmvnorm(c(1, log(0.125), (log(4) / 1.96)^2)),
   prior_is_EXNEX_comp = c(FALSE, FALSE),
   prior_is_EXNEX_inter = FALSE,
   prior_EX_prob_comp = matrix(1,
     nrow = dims$num_groups,
     ncol = 2
   ),
   prior_EX_prob_inter = matrix(1,
     nrow = nlevels(dose_info_combo2$group_id),
     ncol = 1
   ),
   prior_tau_dist = 1
 )
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: 11 out of 42 ewoc metrics have not converged (some Rhats are > 1.1).
#> Be careful when analysing the results! It is recommended to run
#> more iterations and/or setting stronger priors.
#> You may call "summary(trial, summarize='ewoc_check', ...)" for more diagnostic details.
#> Please call "help('blrm_trial', help_type='summary')" for further documentation.
#> Warning: 16 out of 42 ewoc metrics are within the 95% MCMC error of the decision boundary.
#> Be careful when using the imprecise ewoc estimates! It is recommended to run
#> more iterations and review doses close to critical thresholds.
#> You may call "summary(trial, summarize='ewoc_check', ...)" for more diagnostic details.
#> Please call "help('blrm_trial', help_type='summary')" for further documentation.

# print summary of prior specification
prior_summary(combo2_trial_start)
#> Bayesian Logistic Regression Model with EXchangeability-NonEXchangeability
#> 
#> Mixture configuration
#> ---------------------
#> EXNEX components : 0 
#> component
#>    I(log(drug_A/6)) I(log(drug_B/1500)) 
#>                   0                   0 
#> 
#> EXNEX interactions: 0 
#> interaction
#> I(2 * (drug_A/6 * drug_B/1500)/(1 + drug_A/6 * drug_B/1500)) 
#>                                                            0 
#> 
#> Prior probability for exchangeability per group
#>           component
#> group      I(log(drug_A/6)) I(log(drug_B/1500))
#>   trial_A                 1                   1
#>   trial_B                 1                   1
#>   IIT                     1                   1
#>   trial_AB                1                   1
#> 
#>           interaction
#> group      I(2 * (drug_A/6 * drug_B/1500)/(1 + drug_A/6 * drug_B/1500))
#>   trial_A                                                             1
#>   trial_B                                                             1
#>   IIT                                                                 1
#>   trial_AB                                                            1
#> 
#> EXchangable hyperparameter priors
#> ---------------------------------
#> Component parameters
#> Mean mu_log_beta
#>                            prior weight m_intercept m_log_slope s_intercept s_log_slope  rho
#> component           mix                                                                     
#> I(log(drug_A/6))    comp_1          1.0        -1.4         0.0         2.0         1.0  0.0
#> I(log(drug_B/1500)) comp_1          1.0        -1.4         0.0         2.0         1.0  0.0
#> 
#> Heterogeneity tau_log_beta (log-normal)
#>                                      prior weight m_tau_intercept m_tau_log_slope s_tau_intercept s_tau_log_slope   rho
#> stratum   component           mix                                                                                      
#> stratum_1 I(log(drug_A/6))    comp_1         1.00           -1.39           -2.08            0.71            0.71  0.00
#>           I(log(drug_B/1500)) comp_1         1.00           -1.39           -2.08            0.71            0.71  0.00
#> 
#> Correlation LKJ
#> component
#>    I(log(drug_A/6)) I(log(drug_B/1500)) 
#>                   1                   1 
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
#>                            prior weight m_intercept m_log_slope s_intercept s_log_slope  rho
#> component           mix                                                                     
#> I(log(drug_A/6))    comp_1          1.0        -1.4         0.0         2.0         1.0  0.0
#> I(log(drug_B/1500)) comp_1          1.0        -1.4         0.0         2.0         1.0  0.0
#> 
#> Interaction parameters
#> Mean mu_eta
#>       prior   w m[1] s[1]
#> mix                      
#> comp1       1.0  0.0  1.1

# summarize inference at observed dose levels
summary(combo2_trial_start, "data_prediction")
#> # A tibble: 11 × 17
#>    group_id drug_A drug_B num_patients num_toxicities cohort_time stratum_id
#>    <fct>     <dbl>  <dbl>        <dbl>          <dbl>       <dbl> <fct>     
#>  1 trial_A     3      0              3              0           0 all       
#>  2 trial_A     4.5    0              3              0           0 all       
#>  3 trial_A     6      0              6              0           0 all       
#>  4 trial_A     8      0              3              2           0 all       
#>  5 trial_B     0     33.3            3              0           0 all       
#>  6 trial_B     0     50              3              0           0 all       
#>  7 trial_B     0    100              4              0           0 all       
#>  8 trial_B     0    200              9              0           0 all       
#>  9 trial_B     0    400             15              0           0 all       
#> 10 trial_B     0    800             20              2           0 all       
#> 11 trial_B     0   1120             17              4           0 all       
#> # ℹ 10 more variables: dose_id <int>, mean <dbl>, sd <dbl>, `2.5%` <dbl>,
#> #   `50%` <dbl>, `97.5%` <dbl>, prob_underdose <dbl>, prob_target <dbl>,
#> #   prob_overdose <dbl>, ewoc_ok <lgl>

# summarize inference at specified dose levels
summary(combo2_trial_start, "dose_prediction")
#> # A tibble: 42 × 14
#>    group_id drug_A drug_B dose_id stratum_id   mean     sd   `2.5%`   `50%`
#>    <fct>     <dbl>  <dbl>   <int> <fct>       <dbl>  <dbl>    <dbl>   <dbl>
#>  1 trial_A     3        0       1 all        0.0142 0.0138 0.000137 0.00842
#>  2 trial_A     4.5      0       2 all        0.0522 0.0379 0.00545  0.0498 
#>  3 trial_A     6        0       3 all        0.144  0.0903 0.0652   0.0982 
#>  4 trial_A     8        0       4 all        0.373  0.202  0.141    0.380  
#>  5 IIT         0      400       5 all        0.0445 0.0727 0.000623 0.00830
#>  6 IIT         0      600       6 all        0.0641 0.0925 0.00287  0.0183 
#>  7 IIT         0      800       7 all        0.0869 0.107  0.00821  0.0366 
#>  8 IIT         3        0       8 all        0.0212 0.0269 0.000382 0.0136 
#>  9 IIT         3      400       9 all        0.0728 0.112  0.00518  0.0318 
#> 10 IIT         3      600      10 all        0.0967 0.142  0.00701  0.0407 
#> # ℹ 32 more rows
#> # ℹ 5 more variables: `97.5%` <dbl>, prob_underdose <dbl>, prob_target <dbl>,
#> #   prob_overdose <dbl>, ewoc_ok <lgl>


# Update again with new data

# using update() with data argument supplied
# dem <- update(combo2_trial_start, data = codata_combo2)

# alternate way using update() with add_data argument for
# new observations only (those collected after the trial
# design stage).
new_data <- filter(codata_combo2, cohort_time > 0)

combo2_trial <- update(combo2_trial_start, add_data = new_data)
#> stratum_id not given, but only one stratum defined. Assigning first stratum.
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: 17 out of 42 ewoc metrics have not converged (some Rhats are > 1.1).
#> Be careful when analysing the results! It is recommended to run
#> more iterations and/or setting stronger priors.
#> You may call "summary(trial, summarize='ewoc_check', ...)" for more diagnostic details.
#> Please call "help('blrm_trial', help_type='summary')" for further documentation.
#> Warning: 11 out of 42 ewoc metrics are within the 95% MCMC error of the decision boundary.
#> Be careful when using the imprecise ewoc estimates! It is recommended to run
#> more iterations and review doses close to critical thresholds.
#> You may call "summary(trial, summarize='ewoc_check', ...)" for more diagnostic details.
#> Please call "help('blrm_trial', help_type='summary')" for further documentation.

summary(combo2_trial, "data") # cohort_time is tracked
#> # A tibble: 27 × 8
#>    group_id drug_A drug_B num_patients num_toxicities cohort_time stratum_id
#>    <fct>     <dbl>  <dbl>        <dbl>          <dbl>       <dbl> <fct>     
#>  1 trial_A     3      0              3              0           0 all       
#>  2 trial_A     4.5    0              3              0           0 all       
#>  3 trial_A     6      0              6              0           0 all       
#>  4 trial_A     8      0              3              2           0 all       
#>  5 trial_B     0     33.3            3              0           0 all       
#>  6 trial_B     0     50              3              0           0 all       
#>  7 trial_B     0    100              4              0           0 all       
#>  8 trial_B     0    200              9              0           0 all       
#>  9 trial_B     0    400             15              0           0 all       
#> 10 trial_B     0    800             20              2           0 all       
#> # ℹ 17 more rows
#> # ℹ 1 more variable: dose_id <int>
summary(combo2_trial, "data_prediction")
#> # A tibble: 27 × 17
#>    group_id drug_A drug_B num_patients num_toxicities cohort_time stratum_id
#>    <fct>     <dbl>  <dbl>        <dbl>          <dbl>       <dbl> <fct>     
#>  1 trial_A     3      0              3              0           0 all       
#>  2 trial_A     4.5    0              3              0           0 all       
#>  3 trial_A     6      0              6              0           0 all       
#>  4 trial_A     8      0              3              2           0 all       
#>  5 trial_B     0     33.3            3              0           0 all       
#>  6 trial_B     0     50              3              0           0 all       
#>  7 trial_B     0    100              4              0           0 all       
#>  8 trial_B     0    200              9              0           0 all       
#>  9 trial_B     0    400             15              0           0 all       
#> 10 trial_B     0    800             20              2           0 all       
#> # ℹ 17 more rows
#> # ℹ 10 more variables: dose_id <int>, mean <dbl>, sd <dbl>, `2.5%` <dbl>,
#> #   `50%` <dbl>, `97.5%` <dbl>, prob_underdose <dbl>, prob_target <dbl>,
#> #   prob_overdose <dbl>, ewoc_ok <lgl>
summary(combo2_trial, "dose_prediction")
#> # A tibble: 42 × 14
#>    group_id drug_A drug_B dose_id stratum_id     mean       sd    `2.5%`   `50%`
#>    <fct>     <dbl>  <dbl>   <int> <fct>         <dbl>    <dbl>     <dbl>   <dbl>
#>  1 trial_A     3        0       1 all        0.000317 0.000580   6.95e-6 6.01e-5
#>  2 trial_A     4.5      0       2 all        0.00616  0.00688    9.10e-4 3.03e-3
#>  3 trial_A     6        0       3 all        0.0624   0.0316     2.69e-2 6.29e-2
#>  4 trial_A     8        0       4 all        0.482    0.0967     3.14e-1 4.72e-1
#>  5 IIT         0      400       5 all        0.0915   0.0849     7.11e-3 6.60e-2
#>  6 IIT         0      600       6 all        0.168    0.130      3.31e-2 1.29e-1
#>  7 IIT         0      800       7 all        0.258    0.157      8.15e-2 2.39e-1
#>  8 IIT         3        0       8 all        0.000456 0.000550   8.57e-6 2.21e-4
#>  9 IIT         3      400       9 all        0.102    0.0682     1.65e-2 8.36e-2
#> 10 IIT         3      600      10 all        0.213    0.0769     9.70e-2 2.14e-1
#> # ℹ 32 more rows
#> # ℹ 5 more variables: `97.5%` <dbl>, prob_underdose <dbl>, prob_target <dbl>,
#> #   prob_overdose <dbl>, ewoc_ok <lgl>

rm(dims, new_data)

## Recover user set sampling defaults
options(.user_mc_options)
```
