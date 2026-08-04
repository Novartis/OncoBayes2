# Single-agent example

Example data from the application in Neuenschwander, et. al. 2008, from
an "open-label, multicenter, non-comparative, dose-escalation cancer
trial to characterize the safety, tolerability, and pharmacokinetic
profile of a drug and to determine its MTD."

## Usage

``` r
hist_SA
```

## Format

A data frame with 5 rows and 4 variables:

- group_id:

  study

- drug_A:

  dose

- num_patients:

  number of patients

- num_toxicities:

  number of events

## References

Neuenschwander B, Branson M, Gsponer T (2008). “Critical aspects of the
Bayesian approach to phase I cancer trials.” *Statistics in Medicine*,
**27**(13), 2420–2439.
[doi:10.1002/sim.3230](https://doi.org/10.1002/sim.3230) .

## Examples

``` r
.user_mc_options <- options()

## Example from Neuenschwander, B., et al. (2009). Stats in Medicine

dref <- 50

## Since there is no prior information the hierarchical model
## is not used in this example by setting tau to (almost) 0.
blrmfit <- blrm_exnex(
  cbind(num_toxicities, num_patients - num_toxicities) ~
      1 + log(drug_A / dref) |
      0 |
      group_id,
  data = hist_SA,
  prior_EX_mu_comp = mixmvnorm(c(1, logit(1 / 2), log(1), diag(c(2^2, 1)))),
  ## Setting prior_tau_dist=NULL disables the hierarchical prior which is
  ## not required in this example as we analyze a single trial.
  prior_tau_dist = NULL,
  prior_PD = FALSE
)
## Recover user set sampling defaults
options(.user_mc_options)
```
