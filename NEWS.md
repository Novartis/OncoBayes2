# OncoBayes2 0.9-4 - October 17th, 2025

## Enhancements

* Releae of `OncoBayes2` to GitHub along with new `pkgdown`
  [web-site](https://opensource.nibr.com/OncoBayes2/)

# OncoBayes2 0.9-3 - April 25th, 2025

## Bug fixes

* Fix issue with using multiple strata and the new argument
  `prior_tau_dist=NULL` of `blrm_exnex`.
* Fix issues with draws extraction functions which ignored the
  `inc_warmup` argument.

# OncoBayes2 0.9-2 - April 7th, 2025

## Enhancements

* Extend `as_draws_*` functions to `blrm_trial` inputs.

## Bug fixes

* Ensure that the new option `prior_tau_dist=NULL` is forwarded when
  used via `blrm_trial` to `blrm_exnex`.
* Make sure that `rvars` are exported as arrays when the underyling
  parameter is structured.

# OncoBayes2 0.9-1 - March 17th, 2025

## Enhancements

* Adding new vignette on how to derive MAP priors from historical
  data.
* New `sample_map` argument for `blrm_exnex` function allows to sample
  MAP prior for all defined strata of the model. When set to `TRUE`
  the posterior draws contain `map_log_beta` and `map_eta` which
  correspond to respective MAP priors per stratum for new trials.
* Allow argument `prior_tau_dist` of `blrm_exnex` to be set to
  `NULL`. This will disable the hierarchical model structure entirely
  and this simplifies specifying such models without a hierarchical
  structure (data will then be pooled across groups).

## Bug fixes

* Fix issue whenever different number of mixture components are used
  for prios on drug components.
* Ensure that extracted posterior samples as draws objects have
  labeled array dimensions.
* Small efficiency improvement for Stan model.

# OncoBayes2 0.9-0 - February 28th, 2025

## Enhancements

* Introduce mixture prior arguments for `blrm_exnex`. These replace
  the existing arguments for prior specification, which are now
  deprecated and will be removed in the next main version increase of
  the package.
* Add experimental as_draws* functions allowing to extract samples in
  the `posterior` draws format from `blrmfit` objects.
* Improve error messages when Stan fails to sample.

## Bug fixes

* Added several plot improvements and bug fixes. Fixed an issue
  leading to jagged curves in `plot_toxicity_intervals_stacked()`, and
  another leading to wrongly ordered facets in
  `plot_toxicity_intervals()`. Ensured consistent `bayesplot` themeing
  across all plots.
* Fixed sampling of model prior whenever EXNEX is used.
* Simulation-based calibration performance was increased by roughly
  an order of magnitude.
* Fixed issues for BLRMs without interaction term and `cmdstan` 2.36.0
* Results of model runs with the `cmdstanr` backend are now loaded using
  the `brms` package instead of `rstan`. This fixes issues with loading
  results from `cmdstan` 2.36.0
* The intervals reported by the `summary()` function using the `interval_prob`
  argument by default now include the lower end as a closed interval for
  calls where `predictive = FALSE`. For the default underdosing interval,
  this includes the probability of exactly 0 now.
* Fixed an issue relating to an incompatibility of upcoming R 4.5.0 on
  macOS when using C23 with rstan.

# OncoBayes2 0.8-10 - November 4th, 2024

## Enhancements

* Updated Stan model file syntax to use new array syntax as required
  by Stan >=2.33. This upgrades the minimal Stan version to 2.26.
* Start compressing plots in vignette to slim down file size.

# OncoBayes2 0.8-9 - July 20th, 2023

## Bug fixes

* resolve compilation from source issues on some platforms triggered by
  changes in rstantools

# OncoBayes2 0.8-8 - March 3rd, 2023

## Enhancements

* upon package load `OncoBayes2` will now report the date of the release
  and the respective git commit hash used to create the sources of the
  package.

## Bug fixes

* ensure compatibility with posterior 1.4.0 (posterior returns `num()`
  formatted columns, which `OncoBayes2` is not supporting)
* ensure C++17 compatiblity per CRAN (triggers an issue with
  clang 16)

# OncoBayes2 0.8-7 - August 24th, 2022

## Enhancements

* add warning messages when printing `blrmfit` objects for divergent
  transitions or non-convergence of parameters
* add diagnostic extraction functions `nuts_params`, `rhat`, `log_posterior`
  and `neff_ratio` for `blrmfit` objects as defined in `bayesplot`
* speedup `pp_data` by about 25% which is used in all posterior
  summary functions
* add `ewoc_check` summary routine for `blrm_trial` objects. The
  returned summary contains MCMC diagnostics and accuracy estimates of
  the EWOC metric as defined for the trial. This allows to assess the
  MCMC estimation error for EWOC.
* added warning for imprecise estimates of the EWOC metric whenever
  `blrm_trial` objects are printed.

## Bug fixes

* fix broken exclusion of unnecessary random variables whenever
  `save_warmup=FALSE`.

# OncoBayes2 0.8-6 - May 2nd, 2022

* allow cmdstanr as new backend for blrm_exnex. The Stan model file is
  written with the function cmdstanr::write_stan_file to disk such
  that the model binaries can be cached if the user defines the global
  option cmdstanr_write_stan_file_dir. See ?cmdstanr::write_stan_file
  (requires cmdstanr >= 0.5.0)

# OncoBayes2 0.8-5 - February 28th, 2022

* fix issue with example_model not exposing properly objects when
  running with silent=TRUE

# OncoBayes2 0.8-4 - February 18th, 2022

* Switch SBC runs to use clustermq in lieu of batchtools. Also now use
  L'Ecuyer CMG as rng engine during SBC runs.

* add critical_quantile function allowing to calculate critical doses
  which fulfill conditions like EWOC.

# OncoBayes2 0.8-3 - October 7th, 2021

* Make Stan model compile with 2.27+ (lupmf postfix for variables is
  not allowed)

# OncoBayes2 0.8-2 - September 13th, 2021

* Drop unneccessary random variables from the posterior whenever
  save_warmup=FALSE. This decreases the size of the posterior in
  memory by ~60% for the combo2 example.

# OncoBayes2 0.8-1 - September 1st, 2021

* Address CRAN comments

# OncoBayes2 0.8-0 - August 31st, 2021

* Significantly speedup Stan model by dropping normalization of binomial
* The default interaction model for blrm_trial is now saturating
  (blrm_formula_saturating) in logit space, improving behavior for combinations 
  with antagonistic toxicity.
* blrm_trial simplified prior now defaults to ENXEX off, and prior parameter
  defaults have been updated.
* New save_warmup argument to blrm_exnex can be used to disable saving of warmup
  samples, which defaults to the OncoBayes2.MC.save_warmup option, which 
  defaults to TRUE in line with previous versions.
  NOTE: Starting with the next version 0.9.0 the new default will be FALSE!
* Do not crash for summary(trial, "data_prediction") if blrm_trial has NULL or 
  0-row data
* Spurious warning for dose_id NA on blrm_trial instantiation removed
* Spurious warning for dose not prespecified on summary(..., newdata=...)
  removed
* New plot function plot_toxicity_intervals_stacked(), which visualizes model
  posterior or predictive probabilities in discrete intervals, over a continuous
  range of doses
* Automated visual regression testing implemented using the vdiffr package
* Fix in plot_toxicity_curve to avoid overplotting if newdata contains columns
  other than grouping variables and doses (e.g. dose_id)
* Avoid re-parsing model formula whenever newdata is not used in
  respective functions. This prevents that changes to global variables
  (like dref) have any effect whenever the fitted data only is used.

# OncoBayes2 0.7-0 - May 7th, 2021

* Catch data specification errors for wrongly nested groups/strata
* BLRM trials now also print toxicity probability intervals and the EWOC setting
* The blrm_formula_linear and blrm_formula_saturating interaction model formula 
  generators can now be customized to generate specific interactions using the
  specific_interaction_terms argument.
* Add predictive summaries to summary method.
* Switch predictive outputs to be based on averaging sampling density. This avoids
  the need for additional sampling for predictive summaries.
* Add plotting methods for blrmfit and blrm_trial objects.

# OncoBayes2 0.6-5 - May 7th, 2020

* Drop RBesT dependency
* Avoid setting the ggplot2 theme upon package load
* Enhance SBC runs to monitor and report sampler diagnostics

# OncoBayes2 0.6-4 - April 7th, 2020

* Add Operation Qualification script (run-oq.R in inst/extra)
* Add saturating interaction model (blrm_formula_saturating) and associated tests

# OncoBayes2 0.6-3 - March 18th, 2020

* 2nd attempt to work around issues found by clang sanitizer in Stan model.

# OncoBayes2 0.6-2 - March 9th, 2020

* Work around issues found by clang sanitizer in Stan model.

# OncoBayes2 0.6-1 - March 4th, 2020

* Fix issue with upcoming R 4.0 which changes stringsAsFactors default
  to FALSE
* Fix for upcoming tibble 3.0 upgrade which changes conventions of
  growing factors within tibbles

# OncoBayes2 0.6-0 - February 11th, 2020

* Fix issue with model outputs when EXNEX is being used. We recommend
  all users to upgrade to this version. The SBC runs now include
  checks for the per-group estimates to avoid a regression of the bug.
* Correct printing of prior information wrt. to summaries of number of
  strata and groups.

# OncoBayes2 0.5-8 - December 12th, 2019

* Merge data with same group and dosing in blrm_trial before passing
  into blrm_exnex for improved performance
* Ensure consistent sorting of this data to improve reproducibility of
  output with respect to input data permutations
* speedup Stan model which now skips data rows with no cases
* correct blrm_exnex documentation to reflect correctly exchangability
  model used for the interaction model

# OncoBayes2 0.5-7 - December 9th, 2019

* Improve numerical stability of log_inv_logit function, preventing
  NaNs in the output of pp_data and resulting errors in summary()

# OncoBayes2 0.5-6 - November 29th, 2019

* Support summary.blrm_trial(...) syntax passing into summary.blrmfit()

# OncoBayes2 0.5-5 - November 22nd, 2019

* Fixed Roxygen for posterior_predict

# OncoBayes2 0.5-4 - November 18th, 2019

* fix vignette documentation
* correct blrm trial co-data example
* fix test issues per CRAN checks on macosx-old-r run

# OncoBayes2 0.5-3 - November 15th, 2019

* new blrm_trial function which facilitates dose-escalation trial
  conduct by combining key trial design features
* new add_data argument to update function which adds data to existing
  model objects of class blrmfit or blrm_trial
* various smaller bug fixes

# OncoBayes2 0.4-4 - August 28th, 2019

* run in all examples the example code (remove dontrun sections)

# OncoBayes2 0.4-3 - August 27th, 2019

* use message instead of cat in functions using printing except
  summary or print
* suppress by default messages from Stan, can be enable with
  verbose=TRUE
* make examples run with very short sampling

# OncoBayes2 0.4-2 - August 27th, 2019

* Correct loading and exporting of methods
* Add tidybayes example for continuous use and visutalization of model
* Allow multiple central probability widths in prob argument of
  summary method

# OncoBayes2 0.3-0 - July 31st, 2019

* Correct external package loading in examples

# OncoBayes2 0.2-0 - July 27th, 2019

* Added function prior_summary
* Restructured print output
* Added examples single-agent, combo2 and combo3 along with example data sets
* A lot more documentation on reference pages
* New vignette on standard use case of blrm_exnex in Oncology
* Qualified blrm_exnex model with Simulation Based Calibration

# OncoBayes2 0.1-0 - May 15th, 2018

* Initial release
