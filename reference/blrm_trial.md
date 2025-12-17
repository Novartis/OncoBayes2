# Dose-Escalation Trials guided by Bayesian Logistic Regression Model

`blrm_trial` facilitates the conduct of dose escalation studies guided
by Bayesian Logistic Regression Models (BLRM). While the `blrm_exnex`
only fits the BLRM model to data, the `blrm_trial` function standardizes
the specification of the entire trial design and provides various
standardized functions for trial data accrual and derivation of model
summaries needed for dose-escalation decisions.

## Usage

``` r
blrm_trial(
  data,
  dose_info,
  drug_info,
  simplified_prior = FALSE,
  EXNEX_comp = FALSE,
  EX_prob_comp_hist = 1,
  EX_prob_comp_new = 0.8,
  EXNEX_inter = FALSE,
  EX_prob_inter = 1,
  formula_generator = blrm_formula_saturating,
  interval_prob = c(0, 0.16, 0.33, 1),
  interval_max_mass = c(prob_underdose = 1, prob_target = 1, prob_overdose = 0.25),
  ...
)

# S3 method for class 'blrm_trial'
print(x, ...)
```

## Arguments

- data:

  dose-toxicity data available at design stage of trial

- dose_info:

  specificaion of the dose levels as planned for the ongoing trial arms.

- drug_info:

  specification of drugs used in trial arms.

- simplified_prior:

  logical (defaults to `FALSE`) indicating whether a simplified prior
  should be employed based on the `reference_p_dlt` values provided in
  `drug_info`. **Warning:** The simplified prior will change between
  releases. Please read instructions below in the respective section for
  the simplified prior.

- EXNEX_comp:

  logical (default to `TRUE`) indicating whether
  EXchangeable-NonEXchangeable priors should be employed for all
  component parameters

- EX_prob_comp_hist:

  prior weight (\\\[0,1\]\\, default to 1) on exchangeability for the
  component parameters in groups representing historical data

- EX_prob_comp_new:

  prior weight (\\\[0,1\]\\, default to 0.8) on exchangeability for the
  component parameters in groups representing new or concurrent data

- EXNEX_inter:

  logical (default to `FALSE`) indicating whether
  EXchangeable-NonEXchangeable priors should be employed for all
  interaction parameters

- EX_prob_inter:

  prior weight (\\\[0,1\]\\, defaults to 0.8) on exchangeability for the
  interaction parameters

- formula_generator:

  formula generation function (see for example `blrm_formula_linear` or
  `blrm_formula_saturating`). The formula generator defines the employed
  interaction model.

- interval_prob:

  defines the interval probabilities reported in the standard outputs.
  Defaults to `c(0, 0.16, 0.33, 1)`.

- interval_max_mass:

  named vector defining for each interval of the `interval_prob` vector
  a maxmimal admissable probability mass for a given dose level.
  Whenever the posterior probability mass in a given interval exceeds
  the threshold, then the Escalation With Overdose Control (EWOC)
  criterion is considered to be not fullfilled. Dose levels not
  fullfilling EWOC are ineligible for the next cohort of patients. The
  default restricts the overdose probability to less than 0.25.

- ...:

  Additional arguments are forwarded to `blrm_exnex`, i.e. for the
  purpose of prior specification.

- x:

  `blrm_trial` object to print

## Value

The function returns an object of class `blrm_trial`.

## Details

`blrm_trial` constructs an object of class `blrm_trial` which stores the
compelte information about the study design of a dose-escalation trial.
The study design is defined through the data sets (see sections below
for a definition of the columns):

- data (historical data):

  The `data` argument defines available dose-toxicity data at the design
  stage of the trial. Together with the prior of model (without any
  data) this defines the prior used for the trial conduct.

- dose_info:

  Definition of the pre-specified dose levels explored in the ongoing
  trial arms. Thus, all dose-toxcitiy trial data added to the object is
  expected correspond to one of the dose levels in the pre-defined set
  of dose_info.

- drug_info:

  Determines the drugs used in the trial, their units, reference dose
  level and optionally defines the expected probability for a toxicity
  at the reference dose.

Once the `blrm_trial` object is setup the complete trial design is
specified and the model is fitted to the given `data`. This allows
evaluation of the pre-specified dose levels of the trial design wrt. to
safety, i.e. whether the starting dose of the trial fullfills the
escalate with overdose criterion (EWOC) condition.

The `blrm_trial` trial can also be constructed in a 2-step process which
allows for a more convenient specification of the prior since meta data
like number of drugs and the like can be used. See the example section
for details.

After setup of the initial `blrm_trial` object additional data is added
through the use of the `update` method which has a `add_data` argument
intended to add data from the ongoing trial. The `summary` function
finally allows to extract various model summaries. In particular, the
EWOC criterion can be calculated for the pre-defined dose levels of the
trial.

## Methods (by generic)

- `print(blrm_trial)`: print function.

## Simplified prior

As a convenience for the user, a simplified prior can be specifed
whenever the `reference_p_dlt` column is present in the `drug_info` data
set. However, the user is **warned** that the simplified prior will
change in future releases of the package and thus **we strongly
discourage the use of the simplified prior for setting up trial
designs**. The functionality is intended to provide the user a quick
start and as a starting point. The actually instantiated prior can be
seen as demonstrated below in the examples.

## Input data

The data given to the `data` argument of `blrm_trial` is considered as
the available at design stage of the trial. The collected input data
thus does not necessarily need to have the same dose levels as the
pre-specified dose_info for the ongoing trial(s). It's data columns must
include, but are not limited to:

- group_id:

  study

- stratum_id:

  optional, only required for differential discounting of groups

- num_patients:

  number of patients

- num_toxicities:

  number of toxicities

- drug_A:

  Columns for the dose of each treatment component, with column names
  matching the `drug_name` values specified in the `drug_info` argument
  of `blrm_trial`

## Drug info data

The drug information data-set defines drug properties. The fields
included are:

- drug_name:

  name of drug which is also used as column name for the dose

- dose_ref:

  reference dose

- dose_unit:

  units used for drug amounts

- reference_p_dlt:

  optional; if provided, allows setup of a simplified prior

## Dose info data

The `drug_info` data-set pre-specifies the dose levels of the ongoing
trial. Thus, all data added to the `blrm_trial` through the `update`
command must be consistent with the pre-defined dose levels as no other
than those pre-specified ones can be explored in an ongoing trial.

- dose_id:

  optional column which assigns a unique id to each group_id/dose
  combination. If not specified the column is internally generated.

- group_id:

  study

- drug_A:

  Columns for the dose of each treatment component, with column names
  matching the `drug_name` values specified in the `drug_info` argument
  of `blrm_trial`

## References

Babb, J., Rogatko, A., & Zacks, S. (1998). Cancer phase I clinical
trials: efficient dose escalation with overdose control. *Statistics in
medicine*, 17(10), 1103-1120.

Neuenschwander, B., Roychoudhury, S., & Schmidli, H. (2016). On the use
of co-data in clinical trials. *Statistics in Biopharmaceutical
Research*, 8(3), 345-354.

## See also

Other blrm_trial combo2 example:
[`dose_info_combo2`](https://opensource.nibr.com/OncoBayes2/reference/dose_info_combo2.md),
[`drug_info_combo2`](https://opensource.nibr.com/OncoBayes2/reference/drug_info_combo2.md),
[`example-combo2_trial`](https://opensource.nibr.com/OncoBayes2/reference/example-combo2_trial.md)

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
#> Warning: 14 out of 42 ewoc metrics have not converged (some Rhats are > 1.1).
#> Be careful when analysing the results! It is recommended to run
#> more iterations and/or setting stronger priors.
#> You may call "summary(trial, summarize='ewoc_check', ...)" for more diagnostic details.
#> Please call "help('blrm_trial', help_type='summary')" for further documentation.
#> Warning: 19 out of 42 ewoc metrics are within the 95% MCMC error of the decision boundary.
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

# Warning: The simplified prior will change between releases!
# please refer to the combo2_trial example for a complete
# example. You can obtain this example with
# ?example-combo2_trial
# or by running
# example_model("combo2_trial")

## Recover user set sampling defaults
options(.user_mc_options)
```
