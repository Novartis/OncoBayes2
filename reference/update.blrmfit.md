# Update data of a BLRM analysis

Adds data rows to a
[`blrm_exnex()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_exnex.md)
or
[`blrm_trial()`](https://opensource.nibr.com/OncoBayes2/reference/blrm_trial.md)
analysis object.

## Usage

``` r
# S3 method for class 'blrmfit'
update(object, ..., add_data)
```

## Arguments

- object:

  blrmfit analysis object

- ...:

  passed to default `update` command

  The data in `add_data` will be combined with data in `object` using
  `bind_rows`. The indices for groups and stratums (if defined) are
  matched between `add_data` and the data of the analysis `object`.

  Note that the `add_data` argument must be named explicitly as
  demonstrated in the example.

- add_data:

  additional data added to analysis data of `object`

## Examples

``` r
.user_mc_options <- options()

example_model("single_agent", silent = TRUE)

library(tibble)
new_cohort <- tibble(group_id = "trial_A", drug_A = 50, num_patients = 4, num_toxicities = 1)

## this would fail, since add_data argument must be named
## new_blrmfit <- update(blrmfit, new_cohort)
new_blrmfit <- update(blrmfit, add_data = new_cohort)

## Recover user set sampling defaults
options(.user_mc_options)
```
