context("posterior evaluations")

load_combo2_example <- function() {
  load_test_fixture("combo2_example")
}

load_combo3_example <- function() {
  load_test_fixture("combo3_example")
}

load_combo2_trial_example <- function() {
  load_test_fixture("combo2_trial_example")
}

test_that("Outputs of posterior_* functions have expected shapes.", {
  combo2 <- load_combo2_example()
  codata_combo2_alt <- codata_combo2

  num_sim <- nsamples(combo2$blrmfit)

  pp1 <- with(combo2, posterior_predict(blrmfit, newdata = codata_combo2_alt))
  expect_equal(ncol(pp1), nrow(codata_combo2_alt))
  expect_equal(nrow(pp1), num_sim)
  expect_equal(nsamples(combo2$blrmfit), num_sim)
  pp2 <- with(combo2, posterior_linpred(blrmfit, newdata = codata_combo2_alt))
  expect_equal(ncol(pp2), nrow(codata_combo2_alt))
  expect_equal(nrow(pp2), num_sim)
})

test_that("Unkown groups are rejected in posterior_* functions.", {
  combo2 <- load_combo2_example()
  codata_combo2_alt <- codata_combo2
  lev_old <- levels(codata_combo2_alt$group_id)
  levels(codata_combo2_alt$group_id) <- c(
    paste0("new_", lev_old[1]),
    lev_old[-1]
  )

  expect_error(
    with(combo2, posterior_predict(blrmfit, newdata = codata_combo2_alt)),
    regexp = "Found unkown factor levels in grouping: new_trial_A"
  )
  expect_error(
    with(combo2, posterior_linpred(blrmfit, newdata = codata_combo2_alt)),
    regexp = "Found unkown factor levels in grouping: new_trial_A"
  )

  ## same error if the group_id is a character instead
  codata_combo2_alt$group_id <- as.character(codata_combo2_alt$group_id)
  expect_error(
    with(combo2, posterior_predict(blrmfit, newdata = codata_combo2_alt)),
    regexp = "Found unkown factor levels in grouping: new_trial_A"
  )
  expect_error(
    with(combo2, posterior_linpred(blrmfit, newdata = codata_combo2_alt)),
    regexp = "Found unkown factor levels in grouping: new_trial_A"
  )

  ## flip the level definitions
  codata_combo2_alt$group_id <- codata_combo2$group_id
  levels(codata_combo2_alt$group_id)[1:2] <- levels(codata_combo2$group_id)[2:1]
  expect_error(
    with(combo2, posterior_predict(blrmfit, newdata = codata_combo2_alt)),
    regexp = "Mismatch in factor level defintion of grouping"
  )
  expect_error(
    with(combo2, posterior_linpred(blrmfit, newdata = codata_combo2_alt)),
    regexp = "Mismatch in factor level defintion of grouping"
  )
})

test_that("Unkown strata are rejected in posterior_* functions.", {
  combo3 <- load_combo3_example()
  hist_combo3_alt <- hist_combo3
  levels(hist_combo3_alt$stratum_id)[1] <- "BIDflex"

  expect_error(
    with(combo3, posterior_predict(blrmfit, newdata = hist_combo3_alt)),
    regexp = "Found unkown factor levels in stratum: BIDflex"
  )
  expect_error(
    with(combo3, posterior_linpred(blrmfit, newdata = hist_combo3_alt)),
    regexp = "Found unkown factor levels in stratum: BIDflex"
  )

  ## same error if the stratum_id is a character instead
  hist_combo3_alt$stratum_id <- as.character(hist_combo3_alt$stratum_id)
  expect_error(
    with(combo3, posterior_predict(blrmfit, newdata = hist_combo3_alt)),
    regexp = "Found unkown factor levels in stratum: BIDflex"
  )
  expect_error(
    with(combo3, posterior_linpred(blrmfit, newdata = hist_combo3_alt)),
    regexp = "Found unkown factor levels in stratum: BIDflex"
  )

  ## flip the level definitions
  hist_combo3_alt$stratum_id <- hist_combo3$stratum_id
  levels(hist_combo3_alt$stratum_id)[1:2] <- levels(hist_combo3$stratum_id)[2:1]
  expect_error(
    with(combo3, posterior_predict(blrmfit, newdata = hist_combo3_alt)),
    regexp = "Mismatch in factor level defintion of stratum"
  )
  expect_error(
    with(combo3, posterior_linpred(blrmfit, newdata = hist_combo3_alt)),
    regexp = "Mismatch in factor level defintion of stratum"
  )
})


## check blrm_trial objects (basic)
check_trial_posterior_functions <- function(example) {
  with(example, {
    suppressWarnings(suppressMessages(
      trial <- blrm_trial(histdata, dose_info, drug_info)
    ))

    expect_error(posterior_linpred(trial), ".*configure the prior.$")
    expect_error(posterior_interval(trial), ".*configure the prior.$")
    expect_error(predictive_interval(trial), ".*configure the prior.$")
    expect_error(nsamples(trial), ".*configure the prior.$")
  })
}

check_trial_with_prior_posterior_functions <- function(trial) {
  items <- nrow(trial$data)
  draws <- nsamples(trial)
  dims <- c(draws, items)
  expect_equal(dims, dim(posterior_linpred(trial)))
  expect_equal(dims, dim(posterior_predict(trial)))
  expect_equal(posterior_interval(trial$blrmfit), posterior_interval(trial))
  expect_equal(c(items, 2), dim(predictive_interval(trial)))
}

test_that("blrm_trial objects without priors reject posterior methods", {
  check_trial_posterior_functions(examples$single_agent)
  check_trial_posterior_functions(examples$combo2)
})

test_that("blrm_trial objects with priors support posterior methods", {
  single_agent_trial <- load_test_fixture("single_agent_trial")
  combo2_trial <- load_test_fixture("combo2_trial")

  check_trial_with_prior_posterior_functions(single_agent_trial)
  check_trial_with_prior_posterior_functions(combo2_trial)
})


# test S3 methods in alphabetical order
test_that("as_draws and friends have resonable outputs", {
  combo2 <- load_combo2_example()
  draws <- as_draws(
    combo2$blrmfit,
    variable = "mu_log_beta[I(log(drug_A/dref[1])),intercept]"
  )
  expect_s3_class(draws, "draws_list")
  expect_equal(
    variables(draws),
    "mu_log_beta[I(log(drug_A/dref[1])),intercept]"
  )
  expect_equal(ndraws(draws), nsamples(combo2$blrmfit))

  draws <- suppressMessages(as_draws_matrix(
    combo2$blrmfit,
    variable = "mu_log_beta[I(log(drug_A/dref[1])),intercept]"
  ))
  expect_s3_class(draws, "draws_matrix")
  expect_equal(
    variables(draws),
    "mu_log_beta[I(log(drug_A/dref[1])),intercept]"
  )
  expect_equal(ndraws(draws), nsamples(combo2$blrmfit))

  draws <- as_draws_array(
    combo2$blrmfit,
    variable = "mu_log_beta[I(log(drug_A/dref[1])),intercept]"
  )
  expect_s3_class(draws, "draws_array")
  expect_equal(
    variables(draws),
    "mu_log_beta[I(log(drug_A/dref[1])),intercept]"
  )
  expect_equal(ndraws(draws), nsamples(combo2$blrmfit))

  draws <- as_draws_df(
    combo2$blrmfit,
    variable = "mu_log_beta[I(log(drug_A/dref[1])),intercept]"
  )
  expect_s3_class(draws, "draws_df")
  expect_equal(
    variables(draws),
    "mu_log_beta[I(log(drug_A/dref[1])),intercept]"
  )
  expect_equal(ndraws(draws), nsamples(combo2$blrmfit))

  draws <- as_draws_list(
    combo2$blrmfit,
    variable = "mu_log_beta[I(log(drug_A/dref[1])),intercept]"
  )
  expect_s3_class(draws, "draws_list")
  expect_equal(
    variables(draws),
    "mu_log_beta[I(log(drug_A/dref[1])),intercept]"
  )
  expect_equal(ndraws(draws), nsamples(combo2$blrmfit))

  draws <- as_draws_rvars(combo2$blrmfit)
  expect_s3_class(draws, "draws_rvars")
  expect_true(nvariables(draws) > 0)
  expect_equal(ndraws(draws), nsamples(combo2$blrmfit))

  combo2_full <- load_test_fixture("combo2_warmup_blrmfit")
  draws <- as_draws_rvars(combo2_full, inc_warmup = TRUE)
  expect_s3_class(draws, "draws_rvars")
  expect_true(nvariables(draws) > 0)

  ## very_fast_sampling(): 500 iter, 250 warmup, 1 chain
  ## with save_warmup = TRUE: 250 warmup + 250 post-warmup = 500 total
  expect_equal(ndraws(draws), 500L)
  expect_null(combo2_full$stanfit)
})

test_that("as_draws_rvars exports dimension labels", {
  combo2 <- load_combo2_example()
  rv <- as_draws_rvars(combo2$blrmfit, variable = "beta_group")
  ref_dimnames <- list(
    c("trial_A", "trial_B", "IIT", "trial_AB"),
    c("I(log(drug_A/dref[1]))", "I(log(drug_B/dref[2]))"),
    c("intercept", "slope")
  )
  expect_equal(dimnames(rv$beta_group), ref_dimnames)
})

test_that("as_draws_rvars exports expected dimensions for MAP samples", {
  combo_map <- load_test_fixture("combo3_map_blrmfit")
  rv <- as_draws_rvars(combo_map, variable = "map_log_beta")
  expect_equal(dim(rv$map_log_beta), c(2, 3, 2))
})


test_that("as_draws and friends have resonable outputs for blrm_trial", {
  combo2_trial <- load_combo2_trial_example()
  draws <- as_draws(
    combo2_trial$combo2_trial,
    variable = "mu_log_beta[I(log(drug_A/6)),intercept]"
  )
  expect_s3_class(draws, "draws_list")
  expect_equal(variables(draws), "mu_log_beta[I(log(drug_A/6)),intercept]")
  expect_equal(ndraws(draws), nsamples(combo2_trial$combo2_trial))

  draws <- suppressMessages(as_draws_matrix(
    combo2_trial$combo2_trial,
    variable = "mu_log_beta[I(log(drug_A/6)),intercept]"
  ))
  expect_s3_class(draws, "draws_matrix")
  expect_equal(variables(draws), "mu_log_beta[I(log(drug_A/6)),intercept]")
  expect_equal(ndraws(draws), nsamples(combo2_trial$combo2_trial))

  draws <- as_draws_array(
    combo2_trial$combo2_trial,
    variable = "mu_log_beta[I(log(drug_A/6)),intercept]"
  )
  expect_s3_class(draws, "draws_array")
  expect_equal(variables(draws), "mu_log_beta[I(log(drug_A/6)),intercept]")
  expect_equal(ndraws(draws), nsamples(combo2_trial$combo2_trial))

  draws <- as_draws_df(
    combo2_trial$combo2_trial,
    variable = "mu_log_beta[I(log(drug_A/6)),intercept]"
  )
  expect_s3_class(draws, "draws_df")
  expect_equal(variables(draws), "mu_log_beta[I(log(drug_A/6)),intercept]")
  expect_equal(ndraws(draws), nsamples(combo2_trial$combo2_trial))

  draws <- as_draws_list(
    combo2_trial$combo2_trial,
    variable = "mu_log_beta[I(log(drug_A/6)),intercept]"
  )
  expect_s3_class(draws, "draws_list")
  expect_equal(variables(draws), "mu_log_beta[I(log(drug_A/6)),intercept]")
  expect_equal(ndraws(draws), nsamples(combo2_trial$combo2_trial))

  draws <- as_draws_rvars(combo2_trial$combo2_trial)
  expect_s3_class(draws, "draws_rvars")
  expect_true(nvariables(draws) > 0)
  expect_equal(ndraws(draws), nsamples(combo2_trial$combo2_trial))
})

test_that("as_draws_rvars exports dimension labels for blrm_trial", {
  combo2_trial <- load_combo2_trial_example()
  rv <- as_draws_rvars(combo2_trial$combo2_trial, variable = "beta_group")
  ref_dimnames <- list(
    c("trial_A", "trial_B", "IIT", "trial_AB"),
    c("I(log(drug_A/6))", "I(log(drug_B/1500))"),
    c("intercept", "slope")
  )
  expect_equal(dimnames(rv$beta_group), ref_dimnames)
})

test_that("as_draws_rvars exports expected dimensions for MAP samples for blrm_trial", {
  combo_map <- load_test_fixture("combo2_trial_map_blrmfit")
  rv <- as_draws_rvars(combo_map, variable = "map_log_beta")
  expect_equal(dim(rv$map_log_beta), c(1, 2, 2))
})
