context("critical_quantile tests")

## tolerance relates to dose mg scale
eps <- 0.05

suppressPackageStartupMessages(library(dplyr))

test_that("critical interval probabilites are consistent for blrmfit objects", {
  skip_on_cran()

  single_agent_fit <- load_test_fixture("single_agent_example")$blrmfit
  
  ## these tests recover the doses in a data-set by first obtaining
  ## the interval probabilities from the summary method and then
  ## recovering the respective doses

  nd <- mutate(hist_SA, drug_A = drug_A + 5)
  post_inter_1 <- summary(
    single_agent_fit,
    newdata = nd,
    interval_prob = c(0, 0.1, 1)
  )

  ## now recover the doses which the summary corresponds to, lower tails & upper tails
  for (r in seq_len(nrow(nd))) {
    ndr <- nd[r, , drop = FALSE]
    crit_dose_low <- critical_quantile(
      single_agent_fit,
      newdata = ndr,
      x = "drug_A",
      p = post_inter_1[r, "[0,0.1]"],
      qc = 0.1,
      lower.tail = TRUE,
      interval.x = c(0, 20)
    )
    expect_equal(crit_dose_low, ndr$drug_A, tolerance = eps)
    crit_dose_up <- critical_quantile(
      single_agent_fit,
      newdata = ndr,
      x = "drug_A",
      p = post_inter_1[r, "(0.1,1]"],
      qc = 0.1,
      lower.tail = FALSE,
      interval.x = c(0, 20)
    )
    expect_equal(crit_dose_up, ndr$drug_A, tolerance = eps)
  }

  ## check interval probs... these are not entierly unique. Solution
  ## depends on search range, but it's determinstic what you get.
  post_inter_2 <- summary(
    single_agent_fit,
    newdata = nd,
    interval_prob = c(0, 0.1, 0.8, 1)
  )
  for (r in seq_len(nrow(nd))) {
    ndr <- nd[r, , drop = FALSE]
    crit_dose_inter <- critical_quantile(
      single_agent_fit,
      newdata = ndr,
      x = "drug_A",
      p = post_inter_2[r, "(0.1,0.8]"],
      qc = c(0.1, 0.8),
      interval.x = c(2, 30),
      maxiter = 500
    )
    expect_equal(crit_dose_inter, ndr$drug_A, tolerance = eps)
  }
})


test_that("predictive critical interval probabilites are consistent for blrmfit objects", {
  skip_on_cran()

  single_agent_fit <- load_test_fixture("single_agent_example")$blrmfit

  ## these tests recover the doses in a data-set by first obtaining
  ## the interval probabilities from the summary method and then
  ## recovering the respective doses

  nd <- mutate(hist_SA, drug_A = drug_A + 5, num_patients = 10)
  post_inter_1 <- summary(
    single_agent_fit,
    newdata = nd,
    interval_prob = c(-0.1, 0.101, 1),
    predictive = TRUE,
    transform = TRUE
  )

  ## now recover the doses which the summary corresponds to, lower tails & upper tails
  for (r in seq_len(nrow(nd))) {
    ndr <- nd[r, , drop = FALSE]
    crit_dose_low <- critical_quantile(
      single_agent_fit,
      newdata = ndr,
      x = "drug_A",
      p = post_inter_1[r, "(-0.1,0.101]"],
      qc = 0.101,
      lower.tail = TRUE,
      interval.x = c(0, 20),
      predictive = TRUE
    )
    expect_equal(crit_dose_low, ndr$drug_A, tolerance = eps)
    crit_dose_up <- critical_quantile(
      single_agent_fit,
      newdata = ndr,
      x = "drug_A",
      p = post_inter_1[r, "(0.101,1]"],
      qc = 0.10,
      lower.tail = FALSE,
      interval.x = c(0, 20),
      predictive = TRUE
    )
    expect_equal(crit_dose_up, ndr$drug_A, tolerance = eps)
  }

  ## check interval probs... these are not entierly unique. Solution
  ## depends on search range, but it's determinstic what you get.
  post_inter_2 <- summary(
    single_agent_fit,
    newdata = nd,
    interval_prob = c(0, 0.101, 0.79, 1),
    predictive = TRUE,
    transform = TRUE
  )

  for (r in seq_len(nrow(nd))) {
    ndr <- nd[r, , drop = FALSE]
    crit_dose_inter <- critical_quantile(
      single_agent_fit,
      newdata = ndr,
      x = "drug_A",
      p = post_inter_2[r, "(0.101,0.79]"],
      qc = c(0.101, 0.79),
      interval.x = c(2, 30),
      predictive = TRUE
    )
    expect_equal(crit_dose_inter, ndr$drug_A, tolerance = eps)
  }
})


test_that("critical interval probabilites defaults for blrm_trial objects are consistent with standard EWOC", {
  skip_on_cran()

  trial <- load_test_fixture("critical_quantile_combo2_trial")
  dose_info <- select(trial$dose_info, -dose_id)
  dc <- critical_quantile(trial)

  sc <- summary(
    trial,
    newdata = mutate(dose_info, drug1 = dc),
    interval_prob = c(0.33, 1)
  )
  ref <- rep(0.25, times = nrow(dose_info))
  test <- pull(sc, "[0.33,1]")

  expect_equal(test, ref, tolerance = 2 * eps)

  ## trial with non-standard EWOC will trigger an error with defaults
  trial2 <- load_test_fixture("critical_quantile_combo2_trial_nonstandard_ewoc")
  expect_error(critical_quantile(trial2))
})


test_that("critical interval probabilites work for fractionals", {
  skip_on_cran()

  single_agent_fit <- load_test_fixture("single_agent_example")$blrmfit

  ## note:in .model_distribution the labels for the interval_prob
  ## were wrong in this case leading to no output from the
  ## function. Thus, the expectation here is to get output.
  nd <- mutate(hist_SA, drug_A = drug_A + 5)
  crit <- critical_quantile(
    single_agent_fit,
    newdata = nd,
    x = "drug_A",
    p = 0.25,
    qc = c(1 / 3, 1)
  )

  expect_numeric(
    crit,
    lower = 0,
    finite = TRUE,
    any.missing = FALSE,
    all.missing = FALSE,
    len = nrow(nd)
  )
})
