# Summarise model results

Provides model summaries for
[`blrm_exnex()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_exnex.md)
and
[`blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_trial.md)
analyses.

## Usage

``` r
# S3 method for class 'blrmfit'
summary(
  object,
  newdata,
  transform = !predictive,
  prob = 0.95,
  interval_prob,
  predictive = FALSE,
  ...
)
```

## Arguments

- object:

  fitted model object

- newdata:

  optional data frame specifying for what to predict; if missing, then
  the data of the input model `object` is used

- transform:

  logical (defaults to `FALSE`) indicating if the linear predictor on
  the logit link scale is transformed with `inv_logit` to the 0-1
  response scale.

- prob:

  central probability mass to report, i.e. the quantiles 0.5-prob/2 and
  0.5+prob/2 are displayed. Multiple central widths can be specified.

- interval_prob:

  optional vector of sorted quantiles for which the interval
  probabilities are calculated

- predictive:

  logical indicates if the posterior predictive is being summarized.
  Defaults to `FALSE`.

- ...:

  not used in this function

## Value

Returns a `data.frame` of the key summaries of the posterior mean,
standard deviation, central probability interval, median and optional
interval probabilities. Each row of the `data.frame` corresponds to the
respective input data which is by default the same data set as used for
the
[`blrm_exnex()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_exnex.md)
analysis or the data specified in the `newdata` argument.

## Details

The calculated posterior summaries are returned as a `data.frame` and
contain optional interval probabilites for the specified vector of
sorted quantiles. These summaries are calculated on the response scale
by default and can be obtained on the link scale when setting
`transform=FALSE`.

When the results are requested for the predictive distribution with
`predictive=TRUE`, then the link scale refers to the total counts while
the transformed scale divides the (predictive) counts by the number of
trials such that results are on the 0-1 scale.

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

## obtain underdosing (0-0.16], target dosing (0.16-0.33] and
## overdosing (0.33-1] probabilities
summary(blrmfit, interval_prob = c(0, 0.16, 0.33, 1))
#>          mean          sd         2.5%          50%      97.5% [0,0.16]
#> 1 0.002175423 0.004455258 2.102056e-06 5.600326e-05 0.01109878      1.0
#> 2 0.008724355 0.016941761 7.153526e-05 7.739391e-04 0.04382645      1.0
#> 3 0.026667152 0.043294240 1.074384e-03 6.356063e-03 0.11731116      1.0
#> 4 0.100945832 0.099841064 1.588602e-02 5.984106e-02 0.27934022      0.7
#> 5 0.500049422 0.257584238 1.649624e-01 4.445401e-01 0.85698480      0.0
#>   (0.16,0.33] (0.33,1]
#> 1         0.0      0.0
#> 2         0.0      0.0
#> 3         0.0      0.0
#> 4         0.3      0.0
#> 5         0.3      0.7

## obtain predictive distribution for respective cohorts and
## calculate probability for no event, 1 event or >1 event
## note that this does the calculation for the cohort sizes
## as put into the data-set
summary(blrmfit, interval_prob = c(-1, 0, 1, 10), predictive = TRUE)
#>         mean         sd 2.5% 50% 97.5%    (-1,0]       (0,1]       (1,10]
#> 1 0.00652627 0.08135883    0   0     0 0.9935413 0.006391417 6.730507e-05
#> 2 0.03489742 0.19414638    0   0     1 0.9670534 0.031050985 1.895623e-03
#> 3 0.13333576 0.40437497    0   0     1 0.8881358 0.092754269 1.910991e-02
#> 4 0.40378333 0.68606126    0   0     2 0.6942718 0.223424873 8.230338e-02
#> 5 1.00009884 0.78703834    0   1     2 0.3096653 0.380570644 3.097641e-01

## to obtain the predictive for a cohort-size of 6 for all patients
## in the data-set one would need to use the newdata argument, e.g.
summary(blrmfit,
  newdata = transform(hist_SA, num_patients = 6),
  interval_prob = c(-1, 0, 1, 10), predictive = TRUE
)
#>         mean        sd 2.5% 50% 97.5%    (-1,0]      (0,1]       (1,10]
#> 1 0.01305254 0.1164477    0   0     0 0.9872816 0.01238904 0.0003293602
#> 2 0.05234613 0.2442111    0   0     1 0.9523973 0.04312347 0.0044792015
#> 3 0.16000291 0.4542517    0   0     2 0.8710905 0.10227207 0.0266374558
#> 4 0.60567499 0.9020405    0   0     3 0.6070046 0.24112848 0.1518668865
#> 5 3.00029653 1.8142327    0   3     6 0.0889746 0.16092126 0.7501041392

## Recover user set sampling defaults
options(.user_mc_options)
```
