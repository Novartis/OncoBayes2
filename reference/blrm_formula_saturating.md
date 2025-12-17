# Build a BLRM formula with saturating interaction term in logit-space

`blrm_formula_saturating` is a convenience function for generating a
formula for `blrm_trial` and `blrm_exnex` with an interaction of the
form:

\$\$2 \eta \\ \frac{\prod\_{i=1}^N (d_i \big / d_i^\*)}{1 +
\prod\_{i=1}^N (d_i \big / d_i^\*)}\$\$

## Usage

``` r
blrm_formula_saturating(
  ref_doses,
  max_interaction_level = 2,
  specific_interaction_terms = NULL
)
```

## Arguments

- ref_doses:

  Numeric vector of reference doses with names corresponding to drug
  names

- max_interaction_level:

  Highest interaction order to consider \\\[1 - Inf\]\\. Default: 2

- specific_interaction_terms:

  List of custom interaction terms to generate (e.g.
  `list(c("drug1", "drug2"), c("drug1", "drug3"))`).

## Value

The function returns an object of class `blrm_formula`.

## References

Widmer, L.A., Bean, A., Ohlssen, D., Weber, S., Principled Drug-Drug
Interaction Terms for Bayesian Logistic Regression Models of Drug Safety
in Oncology Phase I Combination Trials *arXiv pre-print*, 2023,
[doi:10.48550/arXiv.2302.11437](https://doi.org/10.48550/arXiv.2302.11437)

## Examples

``` r
ref_doses <- c(drug_A = 10, drug_B = 20)

# can be used with blrm_trial
blrm_formula_saturating(ref_doses)
#> $blrm_formula
#> [1] "cbind(num_toxicities, num_patients - num_toxicities) ~ 1 + I(log(drug_A/10)) | 1 + I(log(drug_B/20)) | 0 + I(2 * (drug_A/10 * drug_B/20) / (1 + drug_A/10 * drug_B/20))| stratum_id / group_id"
#> 
#> $num_interaction_terms
#> [1] 1
#> 
#> $num_components
#> [1] 2
#> 
#> $component_names
#> [1] "drug_A" "drug_B"
#> 
#> $max_interaction_level
#> [1] 2
#> 
#> $ref_doses
#> drug_A drug_B 
#>     10     20 
#> 
#> attr(,"class")
#> [1] "blrm_formula"
```
