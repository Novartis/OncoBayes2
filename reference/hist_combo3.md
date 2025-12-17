# Dataset: historical and concurrent data on a three-way combination

This dataset involves a hypothetical dose-escalation study of
combination therapy with three treatment components. From two previous
studies `HistAgent1` and `HistAgent2`, historical data is available on
each of the treatments as single-agents, as well as two of the two-way
combinations. However, due to a difference in treatment schedule between
the `Combo` study and the historical studies, a stratification (through
`stratum_id`) is made between the groups to allow differential
discounting of the alternate-dose data.

## Usage

``` r
hist_combo3
```

## Format

A data frame with 18 rows and 7 variables:

- group_id:

  study

- drug_A:

  dose of Drug A

- drug_B:

  dose of Drug B

- drug_C:

  dose of Drug C

- num_patients:

  number of patients

- num_toxicities:

  number of DLTs

- stratum_id:

  stratum for `group_id`'s used for differential discounting

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)

## example combo3

library(abind)

dref <- c(500, 500, 1000)
num_comp <- 3
num_inter <- choose(3, 2) + 1
num_strata <- nlevels(hist_combo3$stratum_id)
num_groups <- nlevels(hist_combo3$group_id)

blrmfit <- blrm_exnex(
  cbind(num_toxicities, num_patients - num_toxicities) ~
    1 + I(log(drug_A / dref[1])) |
      1 + I(log(drug_B / dref[2])) |
      1 + I(log(drug_C / dref[3])) |
      0
      + I(drug_A / dref[1] * drug_B / dref[2])
        + I(drug_A / dref[1] * drug_C / dref[3])
        + I(drug_B / dref[2] * drug_C / dref[3])
        + I(drug_A / dref[1] * drug_B / dref[2] * drug_C / dref[3]) |
      stratum_id / group_id,
  data = hist_combo3,
  prior_EX_mu_comp = replicate(num_comp, mixmvnorm(c(1, logit(1/3), 0, diag(c(2^2, 1)))), FALSE),
  prior_EX_tau_comp = list(replicate(num_comp,
                                     mixmvnorm(c(1, log(c(0.25, 0.125)),
                                               diag(c(log(4)/1.96, log(4)/1.96)^2))), FALSE),
                           replicate(num_comp,
                                     mixmvnorm(c(1, log(2 * c(0.25, 0.125)),
                                               diag(c(log(4)/1.96, log(4)/1.96)^2))), FALSE)),
  prior_EX_mu_inter = mixmvnorm(c(1, rep.int(0, num_inter),
                                     diag((rep.int(sqrt(2) / 2, num_inter))^2))),
  prior_EX_tau_inter = replicate(num_strata,
                                 mixmvnorm(c(1, rep.int(log(0.25), num_inter),
                                             diag((rep.int(log(2) / 1.96, num_inter))^2))), FALSE),
  prior_EX_prob_comp = matrix(0.9, nrow = num_groups, ncol = num_comp),
  prior_EX_prob_inter = matrix(1.0, nrow = num_groups, ncol = num_inter),
  prior_is_EXNEX_comp = rep(TRUE, num_comp),
  prior_is_EXNEX_inter = rep(FALSE, num_inter),
  prior_tau_dist = 1,
  prior_PD = FALSE
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
## Recover user set sampling defaults
options(.user_mc_options)
```
