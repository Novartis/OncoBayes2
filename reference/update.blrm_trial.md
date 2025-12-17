# Update data and/or prior of a BLRM trial

- Adds data rows to a
  [`blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_trial.md)
  object (add_data argument)

- Replaces data of a
  [`blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_trial.md)
  object (data argument)

- Sets the prior of a
  [`blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_trial.md)
  object (... argument will be passed to blrm_exnex)

## Usage

``` r
# S3 method for class 'blrm_trial'
update(object, ...)
```

## Arguments

- object:

  blrm_trial object

- ...:

  passed to default `update` command of `blrm_exnex`

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)


# the combo2_trial example demonstrates the use of add_data of
# update.blrmfit
example_model("combo2_trial")
#> Running combo2_trial example:
#>  library(tibble)
#>  library(dplyr)
#>  library(tidyr)
#> 
#>  # Combo2 example using blrm_trial functionality
#> 
#>  # construct initial blrm_trial object from built-in example datasets
#>  combo2_trial_setup <- blrm_trial(
#>    data = hist_combo2,
#>    dose_info = dose_info_combo2,
#>    drug_info = drug_info_combo2,
#>    simplified_prior = FALSE
#>  )
#>  
#>  # summary of dimensionality of data structures
#>  dims <- summary(combo2_trial_setup, "dimensionality")
#> 
#>  # Fit the initial model with the historical data and fully specified prior
#> 
#> 
#>  combo2_trial_start <- update(
#>     combo2_trial_setup,
#>     ## bivariate normal prior for drug A and drug B of intercept and
#>     ## log-slope
#>     prior_EX_mu_comp =
#>       replicate(2,
#>                 mixmvnorm(c(1,
#>                             logit(0.2), 0,
#>                             diag(c(2^2, 1))))
#>               , FALSE),
#>     prior_EX_tau_comp =
#>       replicate(2,
#>                 mixmvnorm(c(1,
#>                             log(0.25), log(0.125),
#>                             diag(c(log(4)/1.96, log(4)/1.96)^2)))
#>               , FALSE),
#>     prior_EX_mu_inter = mixmvnorm(c(1, 0, 1.121^2)),
#>     prior_EX_tau_inter = mixmvnorm(c(1, log(0.125), (log(4) / 1.96)^2)),
#>     prior_is_EXNEX_comp = c(FALSE, FALSE),
#>     prior_is_EXNEX_inter = FALSE,
#>     prior_EX_prob_comp = matrix(1,
#>       nrow = dims$num_groups,
#>       ncol = 2
#>     ),
#>     prior_EX_prob_inter = matrix(1,
#>       nrow = nlevels(dose_info_combo2$group_id),
#>       ncol = 1
#>     ),
#>     prior_tau_dist = 1
#>   )
#> 
#>  # print summary of prior specification
#>  prior_summary(combo2_trial_start)
#> 
#>  # summarize inference at observed dose levels
#>  summary(combo2_trial_start, "data_prediction")
#> 
#>  # summarize inference at specified dose levels
#>  summary(combo2_trial_start, "dose_prediction")
#> 
#> 
#>  # Update again with new data
#> 
#>  # using update() with data argument supplied
#>  # dem <- update(combo2_trial_start, data = codata_combo2)
#> 
#>  # alternate way using update() with add_data argument for
#>  # new observations only (those collected after the trial
#>  # design stage).
#>  new_data <- filter(codata_combo2, cohort_time > 0)
#> 
#>  combo2_trial <- update(combo2_trial_start, add_data = new_data)
#> 
#>  summary(combo2_trial, "data") # cohort_time is tracked
#>  summary(combo2_trial, "data_prediction")
#>  summary(combo2_trial, "dose_prediction")
#> 
#>  rm(dims, new_data)
#> 
#> No stratum defined - assigning all groups to single stratum "all"
#> Please configure blrm_exnex using the update() function.
#> Warning: There were 10 transitions after warmup that exceeded the maximum treedepth. Increase max_treedepth above 10. See
#> https://mc-stan.org/misc/warnings.html#maximum-treedepth-exceeded
#> Warning: Examine the pairs() plot to diagnose sampling problems
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
#> Warning: 22 out of 42 ewoc metrics are within the 95% MCMC error of the decision boundary.
#> Be careful when using the imprecise ewoc estimates! It is recommended to run
#> more iterations and review doses close to critical thresholds.
#> You may call "summary(trial, summarize='ewoc_check', ...)" for more diagnostic details.
#> Please call "help('blrm_trial', help_type='summary')" for further documentation.
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
#> Warning: 16 out of 42 ewoc metrics have not converged (some Rhats are > 1.1).
#> Be careful when analysing the results! It is recommended to run
#> more iterations and/or setting stronger priors.
#> You may call "summary(trial, summarize='ewoc_check', ...)" for more diagnostic details.
#> Please call "help('blrm_trial', help_type='summary')" for further documentation.
#> Warning: 7 out of 42 ewoc metrics are within the 95% MCMC error of the decision boundary.
#> Be careful when using the imprecise ewoc estimates! It is recommended to run
#> more iterations and review doses close to critical thresholds.
#> You may call "summary(trial, summarize='ewoc_check', ...)" for more diagnostic details.
#> Please call "help('blrm_trial', help_type='summary')" for further documentation.

## Recover user set sampling defaults
options(.user_mc_options)
```
