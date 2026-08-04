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
.user_mc_options <- options()

example_model("single_agent", silent = TRUE)

## obtain underdosing (0-0.16], target dosing (0.16-0.33] and
## overdosing (0.33-1] probabilities
summary(blrmfit, interval_prob = c(0, 0.16, 0.33, 1))
#>         mean         sd         2.5%          50%      97.5% [0,0.16]
#> 1 0.00701927 0.01829711 3.193827e-08 0.0005916213 0.06277091  0.99875
#> 2 0.01704235 0.03022618 4.196824e-06 0.0046573410 0.10794041  0.99300
#> 3 0.04028952 0.04967560 1.325354e-04 0.0217358570 0.17379939  0.96525
#> 4 0.11463776 0.09407901 3.617333e-03 0.0936283804 0.35028368  0.73875
#> 5 0.44625514 0.20650595 9.260015e-02 0.4399088991 0.83667070  0.08475
#>   (0.16,0.33] (0.33,1]
#> 1     0.00125  0.00000
#> 2     0.00675  0.00025
#> 3     0.03425  0.00050
#> 4     0.22800  0.03325
#> 5     0.24600  0.66925

## obtain predictive distribution for respective cohorts and
## calculate probability for no event, 1 event or >1 event
## note that this does the calculation for the cohort sizes
## as put into the data-set
summary(blrmfit, interval_prob = c(-1, 0, 1, 10), predictive = TRUE)
#>         mean        sd 2.5% 50% 97.5%    (-1,0]      (0,1]      (1,10]
#> 1 0.02105781 0.1513876    0   0     0 0.9800550 0.01887134 0.001073674
#> 2 0.06816939 0.2792281    0   0     1 0.9384946 0.05537610 0.006129270
#> 3 0.20144760 0.4926178    0   0     2 0.8337461 0.13619767 0.030056191
#> 4 0.45855105 0.7156588    0   0     2 0.6528129 0.25470655 0.092480521
#> 5 0.89251028 0.7612431    0   1     2 0.3492674 0.40895489 0.241777697

## to obtain the predictive for a cohort-size of 6 for all patients
## in the data-set one would need to use the newdata argument, e.g.
summary(blrmfit,
  newdata = transform(hist_SA, num_patients = 6),
  interval_prob = c(-1, 0, 1, 10), predictive = TRUE
)
#>         mean        sd 2.5% 50% 97.5%    (-1,0]      (0,1]      (1,10]
#> 1 0.04211562 0.2277301    0   0     1 0.9629396 0.03263697 0.004423426
#> 2 0.10225409 0.3576496    0   0     1 0.9132219 0.07357688 0.013201239
#> 3 0.24173712 0.5531809    0   0     2 0.8087179 0.15016912 0.041112926
#> 4 0.68782657 0.9351123    0   0     3 0.5527840 0.27680859 0.170407412
#> 5 2.67753085 1.6618334    0   3     6 0.1022313 0.17359886 0.724169822

## Recover user set sampling defaults
options(.user_mc_options)
```
