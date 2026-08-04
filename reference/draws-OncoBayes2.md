# Transform `blrmfit` or `blrm_trial` to `draws` objects

Transform a `blrmfit` or `blrm_trial` object to a format supported by
the posterior package.

## Usage

``` r
# S3 method for class 'blrmfit'
as_draws(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'blrmfit'
as_draws_matrix(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'blrmfit'
as_draws_array(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'blrmfit'
as_draws_df(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'blrmfit'
as_draws_list(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'blrmfit'
as_draws_rvars(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'blrm_trial'
as_draws(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'blrm_trial'
as_draws_matrix(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'blrm_trial'
as_draws_array(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'blrm_trial'
as_draws_df(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'blrm_trial'
as_draws_list(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)

# S3 method for class 'blrm_trial'
as_draws_rvars(x, variable = NULL, regex = FALSE, inc_warmup = FALSE, ...)
```

## Arguments

- x:

  A `blrmfit` or `blrm_trial` object.

- variable:

  A character vector providing the variables to extract. By default, all
  variables are extracted.

- regex:

  Logical; Should variable be treated as a (vector of) regular
  expressions? Any variable in `x` matching at least one of the regular
  expressions will be selected. Defaults to `FALSE`.

- inc_warmup:

  Should warmup draws be included? Defaults to `FALSE`.

- ...:

  Arguments passed to individual methods (if applicable).

## Details

To subset iterations, chains, or draws, use the
[`posterior::subset_draws()`](https://mc-stan.org/posterior/reference/subset_draws.html)
method after transforming the input object to a `draws` object.

## See also

[`posterior::draws()`](https://mc-stan.org/posterior/reference/draws.html)
[`posterior::subset_draws()`](https://mc-stan.org/posterior/reference/subset_draws.html)

## Examples

``` r
.user_mc_options <- options()

# fit an example model. See documentation for "combo2" example
example_model("combo2")
#> Running combo2 example:
#>  dref <- c(6, 960)
#> 
#>  num_comp   <- 2 # two investigational drugs
#>  num_inter  <- 1 # one drug-drug interaction needs to be modeled
#>  num_groups <- nlevels(codata_combo2$group_id) # no stratification needed
#>  num_strata <- 1 # no stratification needed
#> 
#>  blrmfit <- blrm_exnex(
#>    cbind(num_toxicities, num_patients - num_toxicities) ~
#>      1 + I(log(drug_A / dref[1])) |
#>        1 + I(log(drug_B / dref[2])) |
#>        0 + I(drug_A / dref[1] * drug_B / dref[2]) |
#>        group_id,
#>    data = codata_combo2,
#>    prior_EX_mu_comp  = list(mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
#>                             mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))),
#>    prior_EX_tau_comp = list(mixmvnorm(c(1,
#>                                         log(0.250), log(0.125),
#>                                         diag(c(log(4)/1.96, log(4)/1.96)^2))),
#>                             mixmvnorm(c(1,
#>                                         log(0.250), log(0.125),
#>                                         diag(c(log(4)/1.96, log(4)/1.96)^2)))),
#>    prior_EX_mu_inter = mixmvnorm(c(1, 0, 1.121^2)),
#>    prior_EX_tau_inter = mixmvnorm(c(1, log(0.125), (log(4) / 1.96)^2)),
#>    prior_is_EXNEX_comp = rep(FALSE, num_comp),
#>    prior_is_EXNEX_inter = rep(FALSE, num_inter),
#>    prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = num_comp),
#>    prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = num_inter),
#>    prior_tau_dist = 1
#>  )

post <- as_draws(blrmfit)

## Recover user set sampling defaults
options(.user_mc_options)
```
