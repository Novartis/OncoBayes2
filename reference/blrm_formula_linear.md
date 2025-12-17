# Build a BLRM formula with linear interaction term in logit-space

`blrm_formula_linear` is a convenience function for generating a formula
for `blrm_trial` and `blrm_exnex` with an interaction of the form:

\$\$\eta \\ \prod\_{i=1}^N (d_i \big / d_i^\*))\$\$

## Usage

``` r
blrm_formula_linear(
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

  List of custom interaction terms to generate (e.g. list(c("drug1",
  "drug2"), c("drug1", "drug3"))). Default: NULL

## Value

The function returns an object of class `blrm_formula`.

## Examples

``` r
ref_doses <- c(drug_A = 10, drug_B = 20)

# can be used with blrm_trial
blrm_formula_linear(ref_doses)
#> $blrm_formula
#> [1] "cbind(num_toxicities, num_patients - num_toxicities) ~ 1 + I(log(drug_A/10)) | 1 + I(log(drug_B/20)) | 0 + I(drug_A/10 * drug_B/20) | stratum_id / group_id"
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
