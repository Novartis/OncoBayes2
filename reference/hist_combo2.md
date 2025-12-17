# Dataset: historical data on two single-agents to inform a combination study

One of two datasets from the application described in Neuenschwander et
al (2016). The risk of DLT is to be studied as a function of dose for
two drugs, drug A and drug B. Historical information on the toxicity
profiles of these two drugs is available from single agent trials
`trial_A` and `trial_B`. A second dataset `codata_combo2` is available
from this application, which includes additional dose-toxicity data from
`trial_AB` and `IIT` of the combination of Drugs A and B.

## Usage

``` r
hist_combo2
```

## Format

A tibble with 11 rows and 5 variables:

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
