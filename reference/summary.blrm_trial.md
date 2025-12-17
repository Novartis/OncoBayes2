# Summarise trial

Provides model summaries for
[`blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_trial.md)
analyses.

## Usage

``` r
# S3 method for class 'blrm_trial'
summary(
  object,
  summarize = c("blrmfit", "blrm_exnex_call", "data", "drug_info", "dose_info",
    "dose_prediction", "data_prediction", "newdata_prediction", "dimensionality",
    "interval_prob", "interval_max_mass", "ewoc_check"),
  ...
)
```

## Arguments

- object:

  [`blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_trial.md)
  object

- summarize:

  one of the following options:

  `blrmfit`

  :   summary of the underlying blrmfit object with further arguments
      ...

  `blrm_exnex_call`

  :   blrm_exnex call used to create the `blrmfit` object

  `drug_info`

  :   drug_info for the trial, contains drugs, reference doses and units

  `dose_info`

  :   dose_info that were defined

  `dose_prediction`

  :   prediction for the defined `dose_info`

  `data`

  :   data that were observed

  `data_prediction`

  :   prediction for the observed data

  `newdata_prediction`

  :   prediction for data provided with the `newdata` argument

  `dimensionality`

  :   numeric vector with entries `"num_components"`,
      `"num_interaction_terms"`, `"num_groups"`, `"num_strata"`

  `interval_prob`

  :   interval probabilities reported in the standard outputs

  `interval_max_mass`

  :   named vector defining for each interval of the `interval_prob`
      vector a maxmimal admissable probability mass for a given dose
      level

  `ewoc_check`

  :   MCMC diagnostic and precision estimates of ewoc defining metrics
      for the defined doses in `dose_info` (default) or for the doses
      provided in the `newdata` argument. Please refer to the details
      for reported diagnostics.

- ...:

  further arguments for
  [`summary.blrmfit()`](https://opensource.nibr.com/OncoBayes2/reference/summary.blrmfit.md)

## Details

The `ewoc_check` summary routine allows to assess the accuracy and
reliability of the ewoc criterion with respect to MCMC sampling noise.
The returned summary provides detailled MCMC convergence and precision
estimates for all criteria defined by `interval_prob` and
`interval_max_mass` which contribute to EWOC metric. That is, for each
interval probability with a maximal mass of less than unity the routine
will return these columns:

- `est`:

  the MCMC estimate defining the critical value. For intervals defined
  by a tail probability this corresponds to the respective critical
  quantile while for interval probabilites this is equal to the interval
  probability.

- `stat`:

  centered and standardized test quantity. The estimate is centered by
  the critical value and scaled by the Monte-Carlo standard error (MCSE)
  of the estimate. Hence, negative (positive) values correspond to the
  constraint being (not) fulfilled. The standardization with the MCSE
  allows to compare the values to standard normal quantiles accordingly.

- `mcse`:

  the Monte-Carlo standard error of the estimate determined with
  [`posterior::mcse_quantile()`](https://mc-stan.org/posterior/reference/mcse_quantile.html)
  (tail probability) or
  [`posterior::mcse_mean()`](https://mc-stan.org/posterior/reference/mcse_mean.html)
  (interval probability) functions.

- `ess`:

  the Monte-Carlo effective sample size of the estimate determined with
  [`posterior::ess_quantile()`](https://mc-stan.org/posterior/reference/ess_quantile.html)
  (tail probability) or
  [`posterior::ess_mean()`](https://mc-stan.org/posterior/reference/ess_mean.html)
  (interval probability) functions.

- `rhat`:

  the Monte-Carlo non-convergence diagnostic Rhat as determined with the
  `[rhat][posterior::rhat] function`.

For the common case of requiring that 33% DLT probability is not
exceeded by more than 25% of the posterior probability mass, the
estimate column `est` contains the 75% quantile \\q\_{75\\}\\ and the
standardized statistic `stat` is defined as:

\$\$\text{stat} = \frac{q\_{75\\} - 33\\}{\text{mcse}\_{q\_{75\\}}}\$\$

The statistic is approximately distributed as a standard normal variate.
The `ewoc_check` summary can be used to ensure that the MCMC estimation
accuracy is sufficient.

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)

# construct initial blrm_trial object from built-in example datasets
combo2_trial_setup <- blrm_trial(
  data = hist_combo2,
  dose_info = dose_info_combo2,
  drug_info = drug_info_combo2,
  simplified_prior = TRUE
)
#> No stratum defined - assigning all groups to single stratum "all"
#> Warning: Simplified prior CAN and WILL change with releases. NOT recommended to use in production. Instantiating a simplified prior - run summary(trial, "blrm_exnex_call") to inspect arguments. 
#> Warning: The largest R-hat is NA, indicating chains have not mixed.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#r-hat
#> Warning: Bulk Effective Samples Size (ESS) is too low, indicating posterior means and medians may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#bulk-ess
#> Warning: Tail Effective Samples Size (ESS) is too low, indicating posterior variances and tail quantiles may be unreliable.
#> Running the chains for more iterations may help. See
#> https://mc-stan.org/misc/warnings.html#tail-ess
#> Warning: 12 out of 42 ewoc metrics have not converged (some Rhats are > 1.1).
#> Be careful when analysing the results! It is recommended to run
#> more iterations and/or setting stronger priors.
#> You may call "summary(trial, summarize='ewoc_check', ...)" for more diagnostic details.
#> Please call "help('blrm_trial', help_type='summary')" for further documentation.
#> Warning: 32 out of 42 ewoc metrics are within the 95% MCMC error of the decision boundary.
#> Be careful when using the imprecise ewoc estimates! It is recommended to run
#> more iterations and review doses close to critical thresholds.
#> You may call "summary(trial, summarize='ewoc_check', ...)" for more diagnostic details.
#> Please call "help('blrm_trial', help_type='summary')" for further documentation.

# extract blrm_call to see setup of the prior as passed to blrm_exnex
summary(combo2_trial_setup, "blrm_exnex_call")
#> blrm_exnex(formula = cbind(num_toxicities, num_patients - num_toxicities) ~ 
#>     1 + I(log(drug_A/6)) | 1 + I(log(drug_B/1500)) | 0 + I(2 * 
#>         (drug_A/6 * drug_B/1500)/(1 + drug_A/6 * drug_B/1500)) | 
#>         stratum_id/group_id, data = data, prior_EX_mu_comp = list(mixmvnorm(c(1, 
#>     logit(0.2), 0, diag(c(2, 0.7)^2))), mixmvnorm(c(1, logit(0.2), 
#>     0, diag(c(2, 0.7)^2)))), prior_EX_tau_comp = replicate(2L, 
#>     mixmvnorm(c(1, log(0.25), log(0.125), diag((c(4, 2)/1.96)^2))), 
#>     FALSE), prior_EX_mu_inter = mixmvnorm(c(1, 0, diag(1.5^2, 
#>     1, 1))), prior_EX_tau_inter = mixmvnorm(c(1, c(log(0.5)), 
#>     diag(c(log(2)/1.96)^2, 1, 1))), prior_is_EXNEX_inter = FALSE, 
#>     prior_is_EXNEX_comp = c(FALSE, FALSE), prior_EX_prob_comp = c(1, 
#>     1, 1, 1, 1, 1, 1, 1), prior_EX_prob_inter = c(1, 1, 1, 1), 
#>     prior_tau_dist = 1)

# extract ewoc precision accuracy
ec <- summary(combo2_trial_setup, "ewoc_check")
#> Warning: 12 out of 42 ewoc metrics have not converged (some Rhats are > 1.1).
#> Be careful when analysing the results! It is recommended to run
#> more iterations and/or setting stronger priors.
#> You may call "summary(trial, summarize='ewoc_check', ...)" for more diagnostic details.
#> Please call "help('blrm_trial', help_type='summary')" for further documentation.
#> Warning: 32 out of 42 ewoc metrics are within the 95% MCMC error of the decision boundary.
#> Be careful when using the imprecise ewoc estimates! It is recommended to run
#> more iterations and review doses close to critical thresholds.
#> You may call "summary(trial, summarize='ewoc_check', ...)" for more diagnostic details.
#> Please call "help('blrm_trial', help_type='summary')" for further documentation.

# find any ewoc metrics which are within 95% MCMC error of the threshold
# these are counted as "imprecise" when printing blrm_trial objects
subset(ec, abs(prob_overdose_stat) < qnorm(0.975))
#> # A tibble: 32 × 10
#>    group_id drug_A drug_B dose_id stratum_id prob_overdose_est
#>    <fct>     <dbl>  <dbl>   <int> <fct>                  <dbl>
#>  1 trial_A     6        0       3 all                    0.272
#>  2 trial_A     8        0       4 all                    0.381
#>  3 IIT         3        0       8 all                    0.133
#>  4 IIT         3      400       9 all                    0.201
#>  5 IIT         3      600      10 all                    0.248
#>  6 IIT         3      800      11 all                    0.311
#>  7 IIT         4.5      0      12 all                    0.173
#>  8 IIT         4.5    400      13 all                    0.256
#>  9 IIT         4.5    600      14 all                    0.310
#> 10 IIT         4.5    800      15 all                    0.366
#> # ℹ 22 more rows
#> # ℹ 4 more variables: prob_overdose_stat <dbl>, prob_overdose_mcse <dbl>,
#> #   prob_overdose_ess <dbl>, prob_overdose_rhat <dbl>

# ensure that the ewoc metric only flags "ok" whenever the MCMC error
# is with 95% below the threshold
ewoc_ok <- ec$prob_overdose_stat < qnorm(0.025)

## Recover user set sampling defaults
options(.user_mc_options)
```
