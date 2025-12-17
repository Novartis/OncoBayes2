# Dataset: historical and concurrent data on a two-way combination

One of two datasets from the application described in Neuenschwander et
al (2016). In the study `trial_AB`, the risk of DLT was studied as a
function of dose for two drugs, drug A and drug B. Historical
information on the toxicity profiles of these two drugs is available
from single agent trials `trial_A` and `trial_B`. Another study `IIT`
was run concurrently to `trial_AB`, and studies the same combination. A
second dataset `hist_combo2` is available from this example, which
includes only the data from the single agent studies, prior to the
initiation of `trial_AB` and `IIT`.

## Usage

``` r
codata_combo2
```

## Format

A data frame with 20 rows and 5 variables:

- group_id:

  study

- drug_A:

  dose of Drug A

- drug_B:

  dose of Drug B

- num_patients:

  number of patients

- num_toxicities:

  number of DLTs

- cohort_time:

  cohort number of patients

## References

Neuenschwander, B., Roychoudhury, S., & Schmidli, H. (2016). On the use
of co-data in clinical trials. *Statistics in Biopharmaceutical
Research*, 8(3), 345-354.

## Examples

``` r
## Setting up dummy sampling for fast execution of example
## Please use 4 chains and 100x more warmup & iter in practice
.user_mc_options <- options(
  OncoBayes2.MC.warmup = 10, OncoBayes2.MC.iter = 20, OncoBayes2.MC.chains = 1,
  OncoBayes2.MC.save_warmup = FALSE
)

dref <- c(6, 960)

num_comp   <- 2 # two investigational drugs
num_inter  <- 1 # one drug-drug interaction needs to be modeled
num_groups <- nlevels(codata_combo2$group_id) # no stratification needed
num_strata <- 1 # no stratification needed

blrmfit <- blrm_exnex(
  cbind(num_toxicities, num_patients - num_toxicities) ~
    1 + I(log(drug_A / dref[1])) |
      1 + I(log(drug_B / dref[2])) |
      0 + I(drug_A / dref[1] * drug_B / dref[2]) |
      group_id,
  data = codata_combo2,
  prior_EX_mu_comp  = list(mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
                           mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))),
  prior_EX_tau_comp = list(mixmvnorm(c(1,
                                       log(0.250), log(0.125),
                                       diag(c(log(4)/1.96, log(4)/1.96)^2))),
                           mixmvnorm(c(1,
                                       log(0.250), log(0.125),
                                       diag(c(log(4)/1.96, log(4)/1.96)^2)))),
  prior_EX_mu_inter = mixmvnorm(c(1, 0, 1.121^2)),
  prior_EX_tau_inter = mixmvnorm(c(1, log(0.125), (log(4) / 1.96)^2)),
  prior_is_EXNEX_comp = rep(FALSE, num_comp),
  prior_is_EXNEX_inter = rep(FALSE, num_inter),
  prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = num_comp),
  prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = num_inter),
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
## Recover user set sampling defaults
options(.user_mc_options)
```
