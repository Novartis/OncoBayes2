# Critical quantile

Calculates the critical value of a model covariate for which a fixed
quantile of the crude rate is equal to the specified (tail or interval)
probability. The specified quantile can relate to the posterior or the
predictive distribution of the response rate.

## Usage

``` r
critical_quantile(object, ...)

# S3 method for class 'blrmfit'
critical_quantile(
  object,
  newdata,
  x,
  p,
  qc,
  lower.tail,
  interval.x,
  extendInt.x = c("auto", "no", "yes", "downX", "upX"),
  log.x,
  predictive = FALSE,
  maxiter = 100,
  ...
)

# S3 method for class 'blrm_trial'
critical_quantile(
  object,
  newdata,
  x,
  p,
  qc,
  lower.tail,
  interval.x,
  extendInt.x = c("auto", "no", "yes", "downX", "upX"),
  log.x,
  predictive = FALSE,
  maxiter = 100,
  ...
)
```

## Arguments

- object:

  fitted model object

- ...:

  not used in this function

- newdata:

  optional data frame specifying for what to predict; if missing, then
  the data of the input model `object` is used

- x:

  character giving the covariate name which is being search for to meet
  the critical conditions. This also supports 'tidy' parameter selection
  by specifying `x = vars(...)`, where `...` is specified the same way
  as in
  [`dplyr::select()`](https://dplyr.tidyverse.org/reference/select.html)
  and similar functions. Examples of using `x` in this way can be found
  in the examples. For `blrm_trial` methods, it defaults to the first
  entry in `summary(blrm_trial, "drug_info")$drug_name`.

- p:

  cumulative probability corresponding to the critical quantile range.

- qc:

  if given as a sorted vector of length 2, then the two entries define
  the interval on the response rate scale which must correspond to the
  cumulative probability `p`. In this case `lower.tail` must be `NULL`
  or left missing. If `qc` is only a single number, then `lower.tail`
  must be set to complete the specification of the tail area.

- lower.tail:

  defines if probabilities are lower or upper tail areas. Must be
  defined if `qc` is a single number and must not be defined if `qc`
  denotes an interval.

- interval.x:

  interval the covariate `x` is searched. Defaults to the minimal and
  maximal value of `x` in the data set used. Whenever `lower.tail` is
  set such that it is known that a tail area is used, then the function
  will automatically attempt to enlarge the search interval since the
  direction of the inverted function is determined around the critical
  value.

- extendInt.x:

  controls if the interval given in `interval.x` should be extended
  during the root finding. The default `"auto"` attempts to guess an
  adequate setting in cases thats possible. Other possible values are
  the same as for the `extendInt` argument of the
  [`stats::uniroot()`](https://rdrr.io/r/stats/uniroot.html) function.

- log.x:

  determines if during the numerical search for the critical value the
  covariate space is logarithmized. By default `x` is log transformed
  whenever all values are positive.

- predictive:

  logical indicates if the critical quantile relates to the posterior
  (`FALSE`) or the predictive (`TRUE`) distribution. Defaults to
  `FALSE`.

- maxiter:

  maximal number of iterations for root finding. Defaults to `100`.

## Value

Vector of length equal to the number of rows of the input data set with
the crticial value for the covariate specified as `x` which fullfills
the requirements as detailled above. May return `NA` for cases where no
solution is found on the specified interval `interval.x`.

## Details

The function searches for a critical value \\x_c\\ of the covariate
\\x\\ (`x`) at which the integral over the model response in the
specified interval \\\pi_1\\ to \\\pi_2\\ (`qc`) is equal to the given
probability \\p\\ (`p`). The given interval is considered a tail area
when `lower.tail` is used such that \\\pi_1 = 0\\ for `TRUE` and
\\\pi_2=1\\ for `FALSE`. At the critical value \\x_c\\ the equality
holds

\$\$p = \int\_{\pi_1}^{\pi_2} p(\pi \| x_c) \\ d\pi .\$\$

Note that a solution not guranteed to exist under all circumstances.
However, for a single agent model and when requesting a tail probability
then a solution does exist. In case an interval probability is
specified, then the solution may not necessarily exist and the root
finding is confined to the range specified in `interval.x`.

In case the solution is requested for the predictive distribution
(`predictive=TRUE`), then the respective problem solved leads to the
equality of

\$\$p = \sum\_{y= r_1 = \lceil n \\ \pi_1 \rceil }^{r_2 = \lceil n \\
\pi_2 \rceil } \int p(y\|\pi, n) \\ p(\pi \| x_c) \\ d\pi .\$\$

Furthermore, the covariate space is log transformed by default whenever
all values of the covariate \\x\\ are positive in the data set. Values
of \\0\\ are shifted into the positive domain by the machine precision
to avoid issues with an ill-defined \\\log(0)\\.

For `blrm_trial` objects the default arguments for `p`, `qc` and
`lower.tail` are set to correspond to the highest interval of
`interval_prob` to be constrained by the respective `interval_max_mass`
(which are defined as part of the `blrm_trial` objects). This will in
the default case correspond to the EWOC metric. These defaults are only
applied if `p`, `qc` and `lower.tail` are missing.

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)

# fit an example model. See documentation for "combo2" example
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

# Find dose of drug_A at which EWOC criterium is just fulfilled
data_trial_ab <- subset(codata_combo2, group_id == "trial_AB")
drug_A_crit <- critical_quantile(blrmfit,
  newdata = data_trial_ab,
  x = "drug_A", p = 0.25, qc = 0.33,
  lower.tail = FALSE
)
data_trial_ab$drug_A <- drug_A_crit
summary(blrmfit, newdata = data_trial_ab, interval_prob = c(0, 0.16, 0.33, 1))
#>        mean         sd       2.5%       50%     97.5% [0,0.16] (0.16,0.33]
#> 1 0.2784784 0.10046646 0.15311191 0.2942422 0.4260879      0.1         0.6
#> 2 0.2575160 0.10161998 0.09526054 0.2553971 0.3886394      0.1         0.6
#> 3 0.2784784 0.10046646 0.15311191 0.2942422 0.4260879      0.1         0.6
#> 4 0.2784784 0.10046646 0.15311191 0.2942422 0.4260879      0.1         0.6
#> 5 0.2575160 0.10161998 0.09526054 0.2553971 0.3886394      0.1         0.6
#> 6 0.2769644 0.08214743 0.18177646 0.2522595 0.4188390      0.0         0.7
#> 7 0.2784784 0.10046646 0.15311191 0.2942422 0.4260879      0.1         0.6
#>   (0.33,1]
#> 1      0.3
#> 2      0.3
#> 3      0.3
#> 4      0.3
#> 5      0.3
#> 6      0.3
#> 7      0.3

## Recover user set sampling defaults
options(.user_mc_options)
```
