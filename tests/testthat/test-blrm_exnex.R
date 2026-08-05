context("blrm_exnex tests")


set.seed(12144)

eps <- 1E-4
eps_low <- 0.02

## Lazy fixture loaders — called inside test_that() so that a missing
## fixture only skips the tests that need it, not the entire file.
load_single_agent <- function() load_test_fixture("single_agent_example")
load_combo2 <- function() load_test_fixture("combo2_example")
load_combo3 <- function() load_test_fixture("combo3_example")

suppressPackageStartupMessages(library(rstan))
suppressPackageStartupMessages(library(dplyr))

## perform some basic checks on the shape of outputs
check_model_basic <- function(fit, envir) {
  skip_on_cran()
  data <- fit$data
  num_obs <- nrow(data)
  ss <- summary(fit, interval_prob = c(0, 0.16, 0.33, 1))

  expect_equal(nrow(ss), num_obs)
  expect_equal(ncol(ss), 8)
  expect_equal(sum(is.na(ss)), 0)
  ss2 <- summary(fit, interval_prob = c(0, 1))
  expect_equal(sum((ss2[, "[0,1]"] - 1) < eps), nrow(data))
  ss3 <- summary(
    fit,
    interval_prob = c(-1, 100),
    predictive = TRUE,
    transform = FALSE
  )
  expect_equal(sum((ss3[, "(-1,100]"] - 1) < eps), nrow(data))

  ## check that the median and 50% interval prob align
  s_med <- ss[, "50%"]
  for (i in seq_along(s_med)) {
    i <- 1
    q_med <- s_med[i]
    s_q <- summary(fit, interval_prob = c(0, q_med, 1))
    expect_true(abs(s_q[i, 6] - 0.5) < eps)
    expect_true(abs(s_q[i, 7] - 0.5) < eps)
  }

  iter <- getOption("OncoBayes2.MC.iter")
  warmup <- getOption("OncoBayes2.MC.warmup")
  chains <- getOption("OncoBayes2.MC.chains")
  num_sim <- chains * (iter - warmup)

  plin <- posterior_linpred(fit)
  expect_equal(nrow(plin), num_sim)
  expect_equal(ncol(plin), nrow(data))
  expect_equal(sum(is.na(ss)), 0)

  envir$fit <- fit
  envir$data1 <- data[1, , drop = FALSE]

  suppressWarnings(suppressMessages(capture.output(
    single_obs_fit <- with(envir, update(fit, data = data1))
  )))
  ss1 <- summary(single_obs_fit)
  expect_equal(nrow(ss1), 1)
  expect_equal(ncol(ss1), 5)
  expect_equal(sum(is.na(ss1)), 0)

  plin1 <- posterior_linpred(single_obs_fit)
  expect_equal(nrow(plin1), num_sim)
  expect_equal(ncol(plin1), 1)
}


test_that("blrm_exnex data handling consistency single-agent", {
  single_agent <- load_single_agent()
  check_model_basic(single_agent$blrmfit, single_agent)
})

test_that("blrm_exnex data handling consistency combo2", {
  combo2 <- load_combo2()
  check_model_basic(combo2$blrmfit, combo2)
})

test_that("blrm_exnex data handling consistency combo3", {
  combo3 <- load_combo3()
  check_model_basic(combo3$blrmfit, combo3)
})

test_that("blrm_exnex data handling consistency single-agent with cmdstanr backend", {
  skip_if_no_cmdstan()
  skip_on_cran()
  withr::with_options(list(OncoBayes2.MC.backend = "cmdstanr"), {
    single_agent_cmdstanr <- run_example("single_agent")
    check_model_basic(single_agent_cmdstanr$blrmfit, single_agent_cmdstanr)
  })
})

test_that("blrm_exnex data handling consistency combo2 with cmdstanr backend", {
  skip_if_no_cmdstan()
  skip_on_cran()
  withr::with_options(list(OncoBayes2.MC.backend = "cmdstanr"), {
    combo2_cmdstanr <- run_example("combo2")
    check_model_basic(combo2_cmdstanr$blrmfit, combo2_cmdstanr)
  })
})

test_that("cmdstanr draws exclude raw and auxiliary variables by default", {
  skip_if_no_cmdstan()
  skip_on_cran()
  withr::with_options(list(OncoBayes2.MC.backend = "cmdstanr"), {
    combo2_cmdstanr <- run_example("combo2")
    variables <- posterior::variables(combo2_cmdstanr$blrmfit$draws)
    excluded <- c(
      "log_beta_raw",
      "eta_raw",
      "tau_log_beta_raw",
      "tau_eta_raw",
      "L_corr_log_beta",
      "L_corr_eta",
      "beta",
      "eta",
      "beta_EX_prob",
      "eta_EX_prob"
    )
    excluded_pattern <- paste0(
      "^(",
      paste(excluded, collapse = "|"),
      ")(\\[|$)"
    )
    expect_false(any(grepl(excluded_pattern, variables)))
    expect_true(any(grepl("^beta_group\\[", variables)))
  })
})

test_that("blrm_exnex data handling consistency combo3 with cmdstanr backend", {
  skip_if_no_cmdstan()
  skip_on_cran()
  withr::with_options(list(OncoBayes2.MC.backend = "cmdstanr"), {
    combo3_cmdstanr <- run_example("combo3")
    check_model_basic(combo3_cmdstanr$blrmfit, combo3_cmdstanr)
  })
})

test_that("interval probabilites are consistent", {
  skip_on_cran()
  withr::local_seed(67894356)

  combo2_sens <- load_combo2()
  sens_low <- hist_combo2 %>%
    mutate(num_toxicities = 0, num_patients = 3 * num_patients)
  combo2_sens$sens_low <- sens_low
  suppressWarnings(suppressMessages(
    outp <- capture.output(
      sens_fit_low <- with(
        combo2_sens,
        update(
          blrmfit,
          iter = 11000,
          warmup = 1000,
          chains = 1,
          init = 0,
          data = sens_low
        )
      ),
      type = "message"
    )
  ))
  s_low <- summary(sens_fit_low, interval_prob = c(0, 0.9, 1))
  expect_equal(sum((s_low[, "[0,0.9]"] - 1) < eps), nrow(sens_low))
  expect_equal(sum((s_low[, "(0.9,1]"]) < eps), nrow(sens_low))

  s_low_L <- summary(
    sens_fit_low,
    interval_prob = c(-100, 4, 100),
    transform = FALSE
  )
  expect_equal(sum((s_low_L[, "[-100,4]"] - 1) < eps), nrow(sens_low))
  expect_equal(sum((s_low_L[, "(4,100]"]) < eps), nrow(sens_low))

  s_low_pp_resp <- summary(
    sens_fit_low,
    interval_prob = c(-1, 0.9, 1),
    transform = TRUE,
    predictive = TRUE
  )
  expect_equal(sum((s_low_pp_resp[, "(-1,0.9]"] - 1) < eps), nrow(sens_low))
  expect_equal(sum((s_low_pp_resp[, "(0.9,1]"]) < eps), nrow(sens_low))

  s_low_pp_count <- summary(
    sens_fit_low,
    interval_prob = c(-1, 10, 100),
    transform = FALSE,
    predictive = TRUE
  )
  expect_equal(sum((s_low_pp_count[, "(-1,10]"] - 1) < eps), nrow(sens_low))
  expect_equal(sum((s_low_pp_count[, "(10,100]"]) < eps), nrow(sens_low))

  sens_high <- hist_combo2 %>%
    mutate(num_patients = 20, num_toxicities = num_patients)
  combo2_sens$sens_high <- sens_high
  outp <- suppressWarnings(suppressMessages(capture.output(
    sens_fit_high <- with(
      combo2_sens,
      update(blrmfit, iter = 11000, warmup = 1000, chains = 1, data = sens_high)
    )
  )))
  s_high <- summary(sens_fit_high, interval_prob = c(0, 0.1, 1))
  expect_equal(sum((s_high[, "(0.1,1]"] - 1) < eps), nrow(sens_low))
  expect_equal(sum((s_high[, "[0,0.1]"]) < eps), nrow(sens_low))

  s_high_L <- summary(
    sens_fit_high,
    interval_prob = c(-100, 0, 100),
    transform = FALSE
  )
  expect_equal(sum((s_high_L[, "[-100,0]"]) < eps), nrow(sens_high))
  expect_equal(sum((s_high_L[, "(0,100]"] - 1) < eps), nrow(sens_high))

  s_high_pp_resp <- summary(
    sens_fit_high,
    interval_prob = c(-1, 0.1, 1),
    transform = TRUE,
    predictive = TRUE
  )
  expect_equal(sum((s_high_pp_resp[, "(-1,0.1]"]) < eps_low), nrow(sens_high))
  expect_equal(
    sum((s_high_pp_resp[, "(0.1,1]"] - 1) < eps_low),
    nrow(sens_high)
  )

  s_high_pp_count <- summary(
    sens_fit_high,
    interval_prob = c(-1, 0, 100),
    transform = FALSE,
    predictive = TRUE
  )
  expect_equal(sum((s_high_pp_count[, "(-1,0]"]) < eps_low), nrow(sens_high))
  expect_equal(
    sum((s_high_pp_count[, "(0,100]"] - 1) < eps_low),
    nrow(sens_high)
  )
})

test_that("predictive interval probabilites are correct", {
  skip_on_cran()
  ## as we use a semi-analytic scheme to get more accurate posterior
  ## predictive summaries, we test here against a simulation based
  ## approach
  single_agent_test <- load_single_agent()
  single_agent_test$test_data <- mutate(
    single_agent_test$blrmfit$data,
    num_patients = 100,
    num_toxicities = c(0, 0, 0, 25, 50)
  )
  fit <- with(
    single_agent_test,
    update(blrmfit, iter = 21000, warmup = 1000, chains = 1, data = test_data)
  )
  ndata <- transform(fit$data, num_patients = 100)
  s_pp <- summary(
    fit,
    newdata = ndata,
    prob = 0.5,
    interval_prob = c(-1, 1),
    predictive = TRUE,
    transform = TRUE
  )
  post <- posterior_predict(fit, newdata = ndata) / 100
  nr <- nrow(s_pp)
  expect_equal(sum(abs(colMeans(post) - s_pp[, "mean"]) < eps_low), nr)
  ## expect_equal(sum( ( apply(post, 2, sd) - s_pp[,"sd"] ) < eps_low ), nr) ## rather unstable estimates...skip
  expect_equal(
    sum(
      abs(apply(post, 2, quantile, probs = 0.5, type = 3) - s_pp[, "50%"]) <
        eps_low
    ),
    nr
  )
  expect_equal(
    sum(
      abs(apply(post, 2, quantile, probs = 0.25, type = 3) - s_pp[, "25%"]) <
        eps_low
    ),
    nr
  )
  expect_equal(
    sum(
      abs(apply(post, 2, quantile, probs = 0.75, type = 3) - s_pp[, "75%"]) <
        eps_low
    ),
    nr
  )
  expect_equal(
    sum(
      abs(apply(post, 2, function(x) mean(x <= 1)) - s_pp[, "(-1,1]"]) < eps_low
    ),
    nr
  )

  num <- fit$data$num_patients
  nr <- length(num)
  test <- sweep(predictive_interval(fit, prob = 0.5), 1, num, "/")
  pp <- sweep(posterior_predict(fit), 2, num, "/")
  ref <- t(apply(pp, 2, quantile, c(0.25, 0.75), type = 3))
  expect_equal(sum(abs(test - ref) < eps_low), 2 * nr)
})

## TODO: Add proper runs which are compared against published results,
## e.g. from the book chapter

test_that("interval probabilites are not NaN", {
  combo2 <- load_combo2()
  ss1 <- summary(combo2$blrmfit, interval_prob = c(0.1, 0.4))
  expect_true(all(!is.na(ss1[, 6])))
  ss2 <- summary(
    combo2$blrmfit,
    transform = FALSE,
    interval_prob = logit(c(0.1, 0.4))
  )
  expect_true(all(!is.na(ss2[, 6])))
})

test_that("correctness of numerical stable log1m_exp_max0", {
  b <- -1E-22
  expect_true(
    abs(OncoBayes2:::log1m_exp_max0(b - 1E-5) - log1p(-exp(b - 1E-5))) < 1E-10
  )
  expect_true(
    abs(OncoBayes2:::log1m_exp_max0(b - 1E-3) - log1p(-exp(b - 1E-3))) < 1E-10
  )
  expect_true(
    abs(OncoBayes2:::log1m_exp_max0(b - 1E-1) - log1p(-exp(b - 1E-1))) < 1E-10
  )
  expect_true(
    abs(OncoBayes2:::log1m_exp_max0(b - 1E-0) - log1p(-exp(b - 1E-0))) < 1E-10
  )
})

test_that("all expected posterior quantiles are returned", {
  combo2 <- load_combo2()
  prob <- 0.95
  f <- function() {
    summary(combo2$blrmfit, prob = prob)
  }
  expect_silent(f())
  ss1 <- f()
  p <- unique(c((1 - prob) / 2, (1 - (1 - prob) / 2), 0.5))
  plab <- as.character(100 * p)
  expect_true(all(
    sapply(plab, function(lab) {
      any(grepl(names(ss1), pattern = lab))
    })
  ))

  prob <- c(0.5, 0.95)
  expect_silent(f())
  ss2 <- f()
  p <- unique(c((1 - prob) / 2, (1 - (1 - prob) / 2), 0.5))
  plab <- as.character(100 * p)
  expect_true(all(
    sapply(plab, function(lab) {
      any(grepl(names(ss2), pattern = lab))
    })
  ))
})


query_fit <- function(fit) {
  suppressWarnings(o <- capture_output(print(fit)))
  prior_summary(fit)
  s1 <- summary(fit)
  s2 <- summary(fit, interval_prob = c(0, 0.16, 0.33, 1.0))
  expect_true(ncol(s1) == ncol(s2) - 3)
  pl <- posterior_linpred(fit)
  pi <- posterior_interval(fit)
  pp <- predictive_interval(fit)
}

test_that("blrmfit methods work for single agent models", {
  single_agent <- load_single_agent()
  ## this test is successfull if methods run without errors as we do
  ## not have running examples given they are too costly.
  query_fit(single_agent$blrmfit)
})

test_that("blrmfit methods work for combo2 models", {
  combo2 <- load_combo2()
  ## this test is successfull if methods run without errors as we do
  ## not have running examples given they are too costly.
  query_fit(combo2$blrmfit)
})

test_that("blrmfit methods work for combo3 models", {
  combo3 <- load_combo3()
  ## this test is successfull if methods run without errors as we do
  ## not have running examples given they are too costly.
  query_fit(combo3$blrmfit)
})


test_that("blrm_exnex accepts single-stratum data sets with general prior definition (deprecated interface)", {
  skip_on_cran()

  withr::local_options(lifecycle_verbosity = "quiet")

  num_comp <- 1 # one investigational drug
  num_inter <- 0 # no drug-drug interactions need to be modeled
  num_groups <- nlevels(hist_SA$group_id) # no stratification needed
  num_strata <- 1 # no stratification needed

  dref <- 50

  hist_SA_alt <- mutate(hist_SA, stratum = factor(1))

  ## in case a single stratum is used in the data, the priors for tau should accept
  suppressWarnings(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 + log(drug_A / dref) | 0 | stratum / group_id,
      data = hist_SA_alt,
      prior_EX_mu_mean_comp = matrix(
        c(
          logit(1 / 2), # mean of intercept on logit scale
          log(1)
        ), # mean of log-slope on logit scale
        nrow = num_comp,
        ncol = 2
      ),
      prior_EX_mu_sd_comp = matrix(
        c(
          2, # sd of intercept
          1
        ), # sd of log-slope
        nrow = num_comp,
        ncol = 2
      ),
      prior_EX_tau_mean_comp = array(0, c(1, num_comp, 2)),
      prior_EX_tau_sd_comp = array(1, c(num_comp, 2)),
      prior_EX_prob_comp = matrix(1, nrow = num_comp, ncol = 1),
      prior_tau_dist = 0,
      prior_PD = FALSE,
      cores = 1,
      iter = 10,
      warmup = 5
    )
  )

  expect_true(nrow(summary(blrmfit)) == nrow(hist_SA_alt))
})

test_that("blrm_exnex accepts single-stratum data sets with general prior definition", {
  num_comp <- 1 # one investigational drug
  num_inter <- 0 # no drug-drug interactions need to be modeled
  num_groups <- nlevels(hist_SA$group_id) # no stratification needed
  num_strata <- 1 # no stratification needed

  dref <- 50

  hist_SA_alt <- mutate(hist_SA, stratum = factor(1))

  suppressWarnings(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 + log(drug_A / dref) | 0 | group_id,
      data = hist_SA,
      prior_EX_mu_comp = mixmvnorm(c(1, logit(1 / 2), log(1), diag(c(2^2, 1)))),
      prior_EX_tau_comp = mixmvnorm(c(1, 0, 0, diag(c(1, 1)))),
      prior_EX_prob_comp = matrix(1, nrow = num_comp, ncol = 1),
      prior_tau_dist = 0,
      prior_PD = FALSE,
      cores = 1,
      iter = 10,
      warmup = 5
    )
  )

  expect_true(nrow(summary(blrmfit)) == nrow(hist_SA_alt))
})

test_that("blrm_exnex summaries do not change depending on global variable definitions of dref (deprecated interface)", {
  skip_on_cran()

  withr::local_options(lifecycle_verbosity = "quiet")

  num_comp <- 1 # one investigational drug
  num_inter <- 0 # no drug-drug interactions need to be modeled
  num_groups <- nlevels(hist_SA$group_id) # no stratification needed
  num_strata <- 1 # no stratification needed

  dref <- 50

  hist_SA_alt <- mutate(hist_SA, stratum = factor(1))

  ## in case a single stratum is used in the data, the priors for tau should accept
  suppressWarnings(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 + log(drug_A / dref) | 0 | stratum / group_id,
      data = hist_SA_alt,
      prior_EX_mu_mean_comp = matrix(
        c(
          logit(1 / 2), # mean of intercept on logit scale
          log(1)
        ), # mean of log-slope on logit scale
        nrow = num_comp,
        ncol = 2
      ),
      prior_EX_mu_sd_comp = matrix(
        c(
          2, # sd of intercept
          1
        ), # sd of log-slope
        nrow = num_comp,
        ncol = 2
      ),
      prior_EX_tau_mean_comp = array(0, c(1, num_comp, 2)),
      prior_EX_tau_sd_comp = array(1, c(num_comp, 2)),
      prior_EX_prob_comp = matrix(1, nrow = num_comp, ncol = 1),
      prior_tau_dist = 0,
      prior_PD = FALSE,
      cores = 1,
      iter = 10,
      warmup = 5
    )
  )

  mean1 <- summary(blrmfit)$mean
  pmean1 <- summary(blrmfit, predictive = TRUE)$mean

  dref <- 1000

  mean2 <- summary(blrmfit)$mean
  pmean2 <- summary(blrmfit, predictive = TRUE)$mean

  expect_equal(mean1, mean2)
  expect_equal(pmean1, pmean2)
})

test_that("blrm_exnex summaries do not change depending on global variable definitions of dref", {
  num_comp <- 1 # one investigational drug
  num_inter <- 0 # no drug-drug interactions need to be modeled
  num_groups <- nlevels(hist_SA$group_id) # no stratification needed
  num_strata <- 1 # no stratification needed

  dref <- 50

  hist_SA_alt <- mutate(hist_SA, stratum = factor(1))

  suppressWarnings(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 + log(drug_A / dref) | 0 | group_id,
      data = hist_SA,
      prior_EX_mu_comp = mixmvnorm(c(1, logit(1 / 2), log(1), diag(c(2^2, 1)))),
      prior_EX_tau_comp = mixmvnorm(c(1, 0, 0, diag(c(1, 1)))),
      prior_EX_prob_comp = matrix(1, nrow = num_comp, ncol = 1),
      prior_tau_dist = 0,
      prior_PD = FALSE,
      cores = 1,
      iter = 10,
      warmup = 5
    )
  )

  mean1 <- summary(blrmfit)$mean
  pmean1 <- summary(blrmfit, predictive = TRUE)$mean

  dref <- 1000

  mean2 <- summary(blrmfit)$mean
  pmean2 <- summary(blrmfit, predictive = TRUE)$mean

  expect_equal(mean1, mean2)
  expect_equal(pmean1, pmean2)
})

test_that("blrm_exnex rejects wrongly nested stratum/group combinations in data sets (deprecated interface)", {
  skip_on_cran()

  withr::local_options(lifecycle_verbosity = "quiet")

  hist_data <- tibble(
    group_id = as.factor(c(
      rep("trial_a", 2),
      rep("trial_b", 3),
      rep("trial_c", 1)
    )),
    stratum_id = as.factor(c(
      rep("reg1", 2),
      rep("reg2", 2),
      "reg3",
      rep("reg1", 1)
    )),
    drug = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 0),
    num_toxicities = c(0, 1, 1, 0, 1, 0),
    num_patients = c(2, 6, 3, 4, 9, 29)
  )

  num_comp <- 1
  num_strata <- nlevels(hist_data$stratum_id)
  num_groups <- nlevels(hist_data$group_id)

  expect_error(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 + I(log(drug)) | 0 | stratum_id / group_id,
      data = hist_data,
      prior_EX_mu_mean_comp = matrix(
        c(logit(0.20), 0), # (E(mu_alpha), E(mu_beta))
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_mu_sd_comp = matrix(
        c(2, 1), # (sd(mu_alpha), sd(mu_beta))
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_tau_mean_comp = abind(
        matrix(log(c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 1 reg1
        matrix(log(2 * c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 2 reg2
        matrix(log(2 * c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 3 reg3
        along = 0
      ),
      prior_EX_tau_sd_comp = abind(
        matrix(
          c(log(4) / 1.96, log(2) / 1.96),
          nrow = num_comp,
          ncol = 2,
          TRUE
        ),
        matrix(
          c(log(4) / 1.96, log(2) / 1.96),
          nrow = num_comp,
          ncol = 2,
          TRUE
        ),
        matrix(
          c(log(4) / 1.96, log(2) / 1.96),
          nrow = num_comp,
          ncol = 2,
          TRUE
        ),
        along = 0
      ),
      prior_is_EXNEX_comp = rep(FALSE, num_comp),
      prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = num_comp),
      prior_tau_dist = 1
    ),
    "^Inconsistent.*"
  )
})


test_that("blrm_exnex rejects wrongly nested stratum/group combinations in data sets", {
  hist_data <- tibble(
    group_id = as.factor(c(
      rep("trial_a", 2),
      rep("trial_b", 3),
      rep("trial_c", 1)
    )),
    stratum_id = as.factor(c(
      rep("reg1", 2),
      rep("reg2", 2),
      "reg3",
      rep("reg1", 1)
    )),
    drug = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 0),
    num_toxicities = c(0, 1, 1, 0, 1, 0),
    num_patients = c(2, 6, 3, 4, 9, 29)
  )

  num_comp <- 1
  num_strata <- nlevels(hist_data$stratum_id)
  num_groups <- nlevels(hist_data$group_id)

  expect_error(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 + log(drug) | 0 | stratum_id / group_id,
      data = hist_data,
      prior_EX_mu_comp = mixmvnorm(c(1, logit(1 / 2), log(1), diag(c(2^2, 1)))),
      prior_EX_tau_comp = replicate(
        3,
        mixmvnorm(c(
          1,
          log(c(0.25, 0.125)),
          diag(c(log(4) / 1.96, log(2) / 1.96)^2)
        )),
        FALSE
      ),
      prior_EX_prob_comp = matrix(1, nrow = num_comp, ncol = 1),
      prior_tau_dist = 1,
      prior_PD = FALSE,
      cores = 1,
      iter = 10,
      warmup = 5
    ),
    "^Inconsistent.*"
  )
})

test_that("update.blrmfit grows the data set", {
  skip_on_cran()

  single_agent <- load_single_agent()
  combo2 <- load_combo2()

  single_agent_new <- single_agent
  single_agent_new$new_cohort_SA <- data.frame(
    group_id = "trial_A",
    num_patients = 4,
    num_toxicities = 2,
    drug_A = 50
  )
  suppressWarnings(
    single_agent_new$new_blrmfit_1 <- with(
      single_agent_new,
      update(blrmfit, add_data = new_cohort_SA)
    )
  )
  expect_true(
    nrow(summary(single_agent_new$new_blrmfit_1)) == nrow(hist_SA) + 1
  )

  ## ensure that the data accumulates
  suppressWarnings(
    new_blrmfit_2 <- with(
      single_agent_new,
      update(new_blrmfit_1, add_data = new_cohort_SA)
    )
  )
  expect_true(nrow(summary(new_blrmfit_2)) == nrow(hist_SA) + 2)

  combo2_new <- combo2
  combo2_new$new_codata <- data.frame(
    group_id = c("IIT", "trial_AB"),
    drug_A = c(8, 8),
    drug_B = c(800, 900),
    drug_C = c(10, 20),
    num_patients = 10,
    num_toxicities = 2,
    stringsAsFactors = TRUE
  )

  ## this one will fail due to a factor levels mismatch
  expect_error(
    new_blrmfit_3 <- with(combo2_new, update(blrmfit, add_data = new_codata)),
    "Mismatch in factor level defintion of grouping",
    fixed = TRUE
  )

  combo2_new$new_codata <- mutate(
    combo2_new$new_codata,
    group_id = as.character(group_id)
  )
  set.seed(123144)
  suppressWarnings(
    new_blrmfit_3 <- with(combo2_new, update(blrmfit, add_data = new_codata))
  )
  expect_true(nrow(summary(new_blrmfit_3)) == nrow(codata_combo2) + 2)
  ## Test that adding dummy data does not change results in the other rows
  set.seed(123144)
  combo2_new_with_dummy <- combo2_new
  combo2_new_with_dummy$new_codata <- add_row(
    combo2_new_with_dummy$new_codata,
    group_id = factor("IIT"),
    drug_A = 1,
    drug_B = 1,
    drug_C = 1,
    num_patients = 0,
    num_toxicities = 0
  )

  suppressWarnings(
    new_blrmfit_3_with_dummy <- with(
      combo2_new_with_dummy,
      update(blrmfit, add_data = new_codata)
    )
  )
  expect_equal(
    nrow(summary(new_blrmfit_3)) + 1,
    nrow(summary(new_blrmfit_3_with_dummy))
  )

  ## test if the log-likelihood is the same for parameter-vector 0
  ## on unconstrained space
  get_stanfit <- function(blrmfit) {
    capture.output(
      model <- sampling(
        OncoBayes2:::stanmodels$blrm_exnex,
        data = blrmfit$standata,
        chains = 0
      )
    )
    model
  }

  new_blrmfit_3_stanfit <- get_stanfit(new_blrmfit_3)
  new_blrmfit_3_with_dummy_stanfit <- get_stanfit(new_blrmfit_3_with_dummy)

  num_pars <- rstan::get_num_upars(new_blrmfit_3_stanfit)
  theta_uconst <- rep(0.0, num_pars)
  log_prob_group <- rstan::log_prob(
    new_blrmfit_3_stanfit,
    theta_uconst,
    gradient = FALSE
  )
  log_prob_group_and_dummy <- rstan::log_prob(
    new_blrmfit_3_with_dummy_stanfit,
    theta_uconst,
    gradient = FALSE
  )
  expect_equal(log_prob_group, log_prob_group_and_dummy)

  ## Same for empty group
  set.seed(123144)
  suppressWarnings(
    new_blrmfit_with_empty_group <- with(
      combo2_new,
      update(blrmfit, data = blrmfit$data)
    )
  )
  set.seed(123144)
  suppressWarnings(
    new_blrmfit_with_empty_group_and_dummy <- with(
      combo2_new,
      update(
        blrmfit,
        data = add_row(
          blrmfit$data,
          group_id = factor("IIT"),
          drug_A = 1,
          drug_B = 1,
          num_patients = 0,
          num_toxicities = 0
        )
      )
    )
  )
  expect_equal(
    nrow(summary(new_blrmfit_with_empty_group)) + 1,
    nrow(summary(new_blrmfit_with_empty_group_and_dummy))
  )
  ## summary(new_blrmfit_with_empty_group) - summary(new_blrmfit_with_empty_group_and_dummy)[1:nrow(summary(new_blrmfit_with_empty_group)),]
  ## change level order and test

  ## test if the log-likelihood is the same for parameter-vector 0
  ## on unconstrained space
  new_blrmfit_with_empty_group_stanfit <- get_stanfit(
    new_blrmfit_with_empty_group
  )
  num_pars <- rstan::get_num_upars(new_blrmfit_with_empty_group_stanfit)
  theta_uconst <- rep(0.0, num_pars)
  combo2_new_blrmfit_stanfit <- get_stanfit(combo2_new$blrmfit)
  log_prob_prob_ref <- rstan::log_prob(
    combo2_new_blrmfit_stanfit,
    theta_uconst,
    gradient = FALSE
  )
  log_prob_empty_group <- rstan::log_prob(
    new_blrmfit_with_empty_group_stanfit,
    theta_uconst,
    gradient = FALSE
  )
  new_blrmfit_with_empty_group_and_dummy_stanfit <- get_stanfit(
    new_blrmfit_with_empty_group_and_dummy
  )
  log_prob_empty_group_and_dummy <- rstan::log_prob(
    new_blrmfit_with_empty_group_and_dummy_stanfit,
    theta_uconst,
    gradient = FALSE
  )
  expect_equal(log_prob_prob_ref, log_prob_empty_group)
  expect_equal(log_prob_prob_ref, log_prob_empty_group_and_dummy)

  combo2_new$flaky_new_codata <- mutate(combo2_new$new_codata, group_id = NULL)

  expect_error(
    new_blrmfit_4 <- with(
      combo2_new,
      update(blrmfit, add_data = flaky_new_codata)
    ),
    "Assertion on 'grouping and/or stratum columns'.*"
  )

  levels_wrong <- levels(codata_combo2$group_id)
  levels_wrong[1:2] <- levels_wrong[2:1]
  combo2_new$wrong_codata <- mutate(
    combo2_new$new_codata,
    group_id = factor(as.character(group_id), levels = levels_wrong)
  )
  expect_true(
    sum(
      levels(combo2_new$wrong_codata$group_id) != levels(codata_combo2$group_id)
    ) ==
      2
  )

  expect_error(
    new_blrmfit_5 <- with(combo2_new, update(blrmfit, add_data = wrong_codata)),
    "Mismatch in factor level defintion of grouping",
    fixed = TRUE
  )
})

test_that("update.blrmfit does regular updating", {
  skip_on_cran()
  single_agent <- load_single_agent()

  single_agent_new <- single_agent
  single_agent_new$only_cohort_SA <- data.frame(
    group_id = "trial_A",
    num_patients = 4,
    num_toxicities = 2,
    drug_A = 50,
    stringsAsFactors = TRUE
  )
  suppressWarnings(
    single_agent_new$only_blrmfit_1 <- with(
      single_agent_new,
      update(blrmfit, data = only_cohort_SA)
    )
  )
  expect_true(nrow(summary(single_agent_new$only_blrmfit_1)) == 1)
})

test_that("update.blrmfit combines data and add_data", {
  skip_on_cran()
  single_agent <- load_single_agent()

  single_agent_new <- single_agent
  single_agent_new$only_cohort_SA <- data.frame(
    group_id = "trial_A",
    num_patients = 4,
    num_toxicities = 2,
    drug_A = 50
  )
  single_agent_new$hist_SA_sub <- hist_SA[1:3, ]
  suppressWarnings(
    single_agent_new$only_blrmfit_1 <- with(
      single_agent_new,
      update(blrmfit, data = hist_SA_sub, add_data = only_cohort_SA)
    )
  )
  expect_true(nrow(summary(single_agent_new$only_blrmfit_1)) == 4)
})

test_that("blrm_exnex properly warns/errors if prior_is_EXNEX is inconsistent from prior_EX_prob (deprecated interface)", {
  skip_on_cran()
  withr::local_options(lifecycle_verbosity = "quiet")

  hist_data <- tibble(
    group_id = as.factor(c(
      rep("trial_a", 2),
      rep("trial_b", 3),
      rep("trial_c", 1)
    )),
    stratum_id = as.factor(c(rep("reg1", 2), rep("reg2", 2), rep("reg2", 2))),
    drug1 = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 0),
    drug2 = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 10),
    num_toxicities = c(0, 1, 1, 0, 1, 0),
    num_patients = c(2, 6, 3, 4, 9, 29)
  )

  num_comp <- 2
  num_strata <- nlevels(hist_data$stratum_id)
  num_groups <- nlevels(hist_data$group_id)
  num_inter <- 1

  expect_error(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug1 / 100)) |
          1 + I(log(drug2 / 100)) |
          0 + I(drug1 / 100 * drug2 / 100) |
          stratum_id / group_id,
      data = hist_data,
      prior_is_EXNEX_comp = rep(FALSE, 2),
      prior_EX_prob_comp = matrix(
        c(0.5, 1, 1),
        nrow = num_groups,
        ncol = num_comp,
        byrow = FALSE
      ), # 0.5 would be ignored
      prior_is_EXNEX_inter = FALSE,
      prior_EX_prob_inter = matrix(
        c(1, 1, 1),
        nrow = num_groups,
        ncol = num_inter
      ),
      prior_EX_mu_mean_inter = rep(0, num_inter),
      prior_EX_mu_sd_inter = rep(1, num_inter),
      prior_EX_tau_mean_inter = matrix(
        log(2) / 1.96,
        nrow = num_strata,
        ncol = num_inter
      ),
      prior_EX_tau_sd_inter = matrix(
        log(2) / 1.96,
        nrow = num_strata,
        ncol = num_inter
      ),
      prior_EX_mu_mean_comp = matrix(
        c(logit(0.20), 0), # (E(mu_alpha), E(mu_beta))
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_mu_sd_comp = matrix(
        c(2, 1), # (sd(mu_alpha), sd(mu_beta))
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_tau_mean_comp = abind(
        matrix(log(c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 1 reg1
        matrix(log(2 * c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 2 reg2
        along = 0
      ),
      prior_EX_tau_sd_comp = abind(
        matrix(
          c(log(4) / 1.96, log(2) / 1.96),
          nrow = num_comp,
          ncol = 2,
          TRUE
        ),
        matrix(
          c(log(4) / 1.96, log(2) / 1.96),
          nrow = num_comp,
          ncol = 2,
          TRUE
        ),
        along = 0
      ),
      prior_tau_dist = 1,
    ),
    "*is_EXNEX*"
  )

  expect_error(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug1 / 100)) |
          1 + I(log(drug2 / 100)) |
          0 + I(drug1 / 100 * drug2 / 100) |
          stratum_id / group_id,
      data = hist_data,
      prior_EX_mu_mean_comp = matrix(
        c(logit(0.20), 0), # (E(mu_alpha), E(mu_beta))
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_mu_sd_comp = matrix(
        c(2, 1), # (sd(mu_alpha), sd(mu_beta))
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_tau_mean_comp = abind(
        matrix(log(c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 1 reg1
        matrix(log(2 * c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 2 reg2
        along = 0
      ),
      prior_EX_tau_sd_comp = abind(
        matrix(
          c(log(4) / 1.96, log(2) / 1.96),
          nrow = num_comp,
          ncol = 2,
          TRUE
        ),
        matrix(
          c(log(4) / 1.96, log(2) / 1.96),
          nrow = num_comp,
          ncol = 2,
          TRUE
        ),
        along = 0
      ),
      prior_is_EXNEX_comp = c(TRUE, FALSE),
      prior_EX_prob_comp = cbind(c(1, 1, 1), c(0.5, 1, 1)), # 0.5 would be ignored
      prior_is_EXNEX_inter = FALSE,
      prior_EX_prob_inter = matrix(
        c(1, 1, 1),
        nrow = num_groups,
        ncol = num_inter
      ),
      prior_EX_mu_mean_inter = rep(0, num_inter),
      prior_EX_mu_sd_inter = rep(1, num_inter),
      prior_EX_tau_mean_inter = matrix(
        log(2) / 1.96,
        nrow = num_strata,
        ncol = num_inter
      ),
      prior_EX_tau_sd_inter = matrix(
        log(2) / 1.96,
        nrow = num_strata,
        ncol = num_inter
      ),
      prior_tau_dist = 1
    ),
    "*is_EXNEX*"
  )

  expect_error(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug1 / 100)) |
          1 + I(log(drug2 / 100)) |
          0 + I(drug1 / 100 * drug2 / 100) |
          stratum_id / group_id,
      data = hist_data,
      prior_EX_mu_mean_comp = matrix(
        c(logit(0.20), 0), # (E(mu_alpha), E(mu_beta))
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_mu_sd_comp = matrix(
        c(2, 1), # (sd(mu_alpha), sd(mu_beta))
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_tau_mean_comp = abind(
        matrix(log(c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 1 reg1
        matrix(log(2 * c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 2 reg2
        along = 0
      ),
      prior_EX_tau_sd_comp = abind(
        matrix(
          c(log(4) / 1.96, log(2) / 1.96),
          nrow = num_comp,
          ncol = 2,
          TRUE
        ),
        matrix(
          c(log(4) / 1.96, log(2) / 1.96),
          nrow = num_comp,
          ncol = 2,
          TRUE
        ),
        along = 0
      ),
      prior_is_EXNEX_comp = c(FALSE, FALSE),
      prior_is_EXNEX_inter = FALSE,
      prior_EX_prob_comp = cbind(c(1, 1, 1), c(1, 1, 1)),
      prior_EX_prob_inter = matrix(
        c(0.5, 1, 1),
        nrow = num_groups,
        ncol = num_inter
      ),
      prior_tau_dist = 1,
      prior_EX_mu_mean_inter = rep(0, num_inter),
      prior_EX_mu_sd_inter = rep(1, num_inter),
      prior_EX_tau_mean_inter = matrix(
        log(2) / 1.96,
        nrow = num_strata,
        ncol = num_inter
      ),
      prior_EX_tau_sd_inter = matrix(
        log(2) / 1.96,
        nrow = num_strata,
        ncol = num_inter
      ),
    ),
    "*is_EXNEX*"
  )

  expect_warning(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug1 / 100)) |
          1 + I(log(drug2 / 100)) |
          0 + I(drug1 / 100 * drug2 / 100) |
          stratum_id / group_id,
      data = hist_data,
      prior_EX_mu_mean_comp = matrix(
        c(logit(0.20), 0), # (E(mu_alpha), E(mu_beta))
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_mu_sd_comp = matrix(
        c(2, 1), # (sd(mu_alpha), sd(mu_beta))
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_tau_mean_comp = abind(
        matrix(log(c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 1 reg1
        matrix(log(2 * c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 2 reg2
        along = 0
      ),
      prior_EX_tau_sd_comp = abind(
        matrix(
          c(log(4) / 1.96, log(2) / 1.96),
          nrow = num_comp,
          ncol = 2,
          TRUE
        ),
        matrix(
          c(log(4) / 1.96, log(2) / 1.96),
          nrow = num_comp,
          ncol = 2,
          TRUE
        ),
        along = 0
      ),
      prior_is_EXNEX_comp = c(TRUE, TRUE),
      prior_EX_prob_comp = cbind(c(1, 1, 1), c(1, 1, 1)), # 0.5 would be ignored
      prior_is_EXNEX_inter = FALSE,
      prior_EX_prob_inter = matrix(
        c(1, 1, 1),
        nrow = num_groups,
        ncol = num_inter
      ),
      prior_EX_mu_mean_inter = rep(0, num_inter),
      prior_EX_mu_sd_inter = rep(1, num_inter),
      prior_EX_tau_mean_inter = matrix(
        log(2) / 1.96,
        nrow = num_strata,
        ncol = num_inter
      ),
      prior_EX_tau_sd_inter = matrix(
        log(2) / 1.96,
        nrow = num_strata,
        ncol = num_inter
      ),
      prior_tau_dist = 1,
      init = 0
    ),
    "*is_EXNEX*"
  )

  expect_warning(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug1 / 100)) |
          1 + I(log(drug2 / 100)) |
          0 + I(drug1 / 100 * drug2 / 100) |
          stratum_id / group_id,
      data = hist_data,
      prior_EX_mu_mean_comp = matrix(
        c(logit(0.20), 0), # (E(mu_alpha), E(mu_beta))
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_mu_sd_comp = matrix(
        c(2, 1), # (sd(mu_alpha), sd(mu_beta))
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_tau_mean_comp = abind(
        matrix(log(c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 1 reg1
        matrix(log(2 * c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 2 reg2
        along = 0
      ),
      prior_EX_tau_sd_comp = abind(
        matrix(
          c(log(4) / 1.96, log(2) / 1.96),
          nrow = num_comp,
          ncol = 2,
          TRUE
        ),
        matrix(
          c(log(4) / 1.96, log(2) / 1.96),
          nrow = num_comp,
          ncol = 2,
          TRUE
        ),
        along = 0
      ),
      prior_is_EXNEX_comp = c(FALSE, FALSE),
      prior_is_EXNEX_inter = TRUE,
      prior_EX_prob_comp = cbind(c(1, 1, 1), c(1, 1, 1)),
      prior_EX_prob_inter = matrix(
        c(1, 1, 1),
        nrow = num_groups,
        ncol = num_inter
      ),
      prior_EX_mu_mean_inter = rep(0, num_inter),
      prior_EX_mu_sd_inter = rep(1, num_inter),
      prior_EX_tau_mean_inter = matrix(
        log(2) / 1.96,
        nrow = num_strata,
        ncol = num_inter
      ),
      prior_EX_tau_sd_inter = matrix(
        log(2) / 1.96,
        nrow = num_strata,
        ncol = num_inter
      ),
      prior_tau_dist = 1,
      init = 0
    ),
    "*is_EXNEX*"
  )
})

test_that("blrm_exnex properly warns/errors if prior_is_EXNEX is inconsistent from prior_EX_prob", {
  skip_on_cran()

  hist_data <- tibble(
    group_id = as.factor(c(
      rep("trial_a", 2),
      rep("trial_b", 3),
      rep("trial_c", 1)
    )),
    stratum_id = as.factor(c(rep("reg1", 2), rep("reg2", 2), rep("reg2", 2))),
    drug1 = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 0),
    drug2 = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 10),
    num_toxicities = c(0, 1, 1, 0, 1, 0),
    num_patients = c(2, 6, 3, 4, 9, 29)
  )

  num_comp <- 2
  num_strata <- nlevels(hist_data$stratum_id)
  num_groups <- nlevels(hist_data$group_id)
  num_inter <- 1

  expect_error(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug1 / 100)) |
          1 + I(log(drug2 / 100)) |
          0 + I(drug1 / 100 * drug2 / 100) |
          group_id,
      data = hist_data,
      prior_EX_mu_comp = list(
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))
      ),
      prior_EX_tau_comp = list(
        mixmvnorm(c(
          1,
          log(0.250),
          log(0.125),
          diag(c(log(4) / 1.96, log(4) / 1.96)^2)
        )),
        mixmvnorm(c(
          1,
          log(2 * 0.250),
          log(2 * 0.125),
          diag(c(log(4) / 1.96, log(4) / 1.96)^2)
        ))
      ),
      prior_EX_mu_inter = mixmvnorm(c(1, 0, 1)),
      prior_EX_tau_inter = mixmvnorm(c(1, log(2) / 1.96, (log(2) / 1.96)^2)),
      prior_is_EXNEX_comp = rep(FALSE, num_comp),
      prior_is_EXNEX_inter = rep(FALSE, num_inter),
      prior_EX_prob_comp = matrix(
        c(0.5, 1, 1),
        nrow = num_groups,
        ncol = num_comp,
        byrow = FALSE
      ), # 0.5 would be ignored
      prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = num_inter),
      prior_tau_dist = 1
    ),
    "*is_EXNEX*"
  )

  expect_error(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug1 / 100)) |
          1 + I(log(drug2 / 100)) |
          0 + I(drug1 / 100 * drug2 / 100) |
          group_id,
      data = hist_data,
      prior_EX_mu_comp = list(
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))
      ),
      prior_EX_tau_comp = list(
        mixmvnorm(c(
          1,
          log(0.250),
          log(0.125),
          diag(c(log(4) / 1.96, log(4) / 1.96)^2)
        )),
        mixmvnorm(c(
          1,
          log(2 * 0.250),
          log(2 * 0.125),
          diag(c(log(4) / 1.96, log(4) / 1.96)^2)
        ))
      ),
      prior_EX_mu_inter = mixmvnorm(c(1, 0, 1)),
      prior_EX_tau_inter = mixmvnorm(c(1, log(2) / 1.96, (log(2) / 1.96)^2)),
      prior_is_EXNEX_comp = rep(FALSE, num_comp),
      prior_is_EXNEX_inter = rep(FALSE, num_inter),
      prior_EX_prob_comp = cbind(c(1, 1, 1), c(0.5, 1, 1)), # 0.5 would be ignored
      prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = num_inter),
      prior_tau_dist = 1
    ),
    "*is_EXNEX*"
  )

  expect_error(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug1 / 100)) |
          1 + I(log(drug2 / 100)) |
          0 + I(drug1 / 100 * drug2 / 100) |
          group_id,
      data = hist_data,
      prior_EX_mu_comp = list(
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))
      ),
      prior_EX_tau_comp = list(
        mixmvnorm(c(
          1,
          log(0.250),
          log(0.125),
          diag(c(log(4) / 1.96, log(4) / 1.96)^2)
        )),
        mixmvnorm(c(
          1,
          log(2 * 0.250),
          log(2 * 0.125),
          diag(c(log(4) / 1.96, log(4) / 1.96)^2)
        ))
      ),
      prior_EX_mu_inter = mixmvnorm(c(1, 0, 1)),
      prior_EX_tau_inter = mixmvnorm(c(1, log(2) / 1.96, (log(2) / 1.96)^2)),
      prior_is_EXNEX_comp = rep(FALSE, num_comp),
      prior_is_EXNEX_inter = rep(FALSE, num_inter),
      prior_EX_prob_comp = cbind(c(1, 1, 1), c(1, 1, 1)),
      prior_EX_prob_inter = matrix(
        c(0.5, 1, 1),
        nrow = num_groups,
        ncol = num_inter
      ), # 0.5 would be ignored
      prior_tau_dist = 1
    ),
    "*is_EXNEX*"
  )

  ## warning on efficiency loss due to enabling EXNEX for the components, but not using it
  expect_warning(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug1 / 100)) |
          1 + I(log(drug2 / 100)) |
          0 + I(drug1 / 100 * drug2 / 100) |
          group_id,
      data = hist_data,
      prior_EX_mu_comp = list(
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))
      ),
      prior_EX_tau_comp = list(
        mixmvnorm(c(
          1,
          log(0.250),
          log(0.125),
          diag(c(log(4) / 1.96, log(4) / 1.96)^2)
        )),
        mixmvnorm(c(
          1,
          log(2 * 0.250),
          log(2 * 0.125),
          diag(c(log(4) / 1.96, log(4) / 1.96)^2)
        ))
      ),
      prior_EX_mu_inter = mixmvnorm(c(1, 0, 1)),
      prior_EX_tau_inter = mixmvnorm(c(1, log(2) / 1.96, (log(2) / 1.96)^2)),
      prior_is_EXNEX_comp = rep(TRUE, num_comp),
      prior_is_EXNEX_inter = rep(FALSE, num_inter),
      prior_EX_prob_comp = cbind(c(1, 1, 1), c(1, 1, 1)),
      prior_EX_prob_inter = matrix(
        c(1, 1, 1),
        nrow = num_groups,
        ncol = num_inter
      ),
      prior_tau_dist = 1,
      init = 0
    ),
    "*is_EXNEX*"
  )

  ## same for the interactions
  expect_warning(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug1 / 100)) |
          1 + I(log(drug2 / 100)) |
          0 + I(drug1 / 100 * drug2 / 100) |
          group_id,
      data = hist_data,
      prior_EX_mu_comp = list(
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))
      ),
      prior_EX_tau_comp = list(
        mixmvnorm(c(
          1,
          log(0.250),
          log(0.125),
          diag(c(log(4) / 1.96, log(4) / 1.96)^2)
        )),
        mixmvnorm(c(
          1,
          log(2 * 0.250),
          log(2 * 0.125),
          diag(c(log(4) / 1.96, log(4) / 1.96)^2)
        ))
      ),
      prior_EX_mu_inter = mixmvnorm(c(1, 0, 1)),
      prior_EX_tau_inter = mixmvnorm(c(1, log(2) / 1.96, (log(2) / 1.96)^2)),
      prior_is_EXNEX_comp = rep(FALSE, num_comp),
      prior_is_EXNEX_inter = rep(TRUE, num_inter),
      prior_EX_prob_comp = cbind(c(1, 1, 1), c(1, 1, 1)),
      prior_EX_prob_inter = matrix(
        c(1, 1, 1),
        nrow = num_groups,
        ncol = num_inter
      ),
      prior_tau_dist = 1,
      init = 0
    ),
    "*is_EXNEX*"
  )
})


test_that("blrm_exnex posterior predictons are not randomly shuffled in their order", {
  single_agent <- load_single_agent()
  ## as some MCMC diagnostics rely on the original order of the MCMC
  ## chain, we test here for the intercept that no permutation is
  ## being done
  post_rv <- as_draws_rvars(single_agent$blrmfit, variable = "beta_group")
  post_pl <- posterior_linpred(
    single_agent$blrmfit,
    newdata = mutate(hist_SA, drug_A = single_agent$dref)[1, ]
  )
  apost_rv <- posterior::draws_of(post_rv$beta_group, TRUE)
  iter <- posterior::niterations(post_rv)
  for (i in seq_len(posterior::nchains(post_rv))) {
    expect_true(all(
      rank(apost_rv[, i, 1, 1, 1]) ==
        rank(post_pl[seq(iter * (i - 1), iter * i - 1) + 1, ])
    ))
    expect_true(
      all(
        abs(
          apost_rv[, i, 1, 1, 1] -
            post_pl[seq(iter * (i - 1), iter * i - 1) + 1, ]
        ) <
          eps
      ),
      "Draws per chain must match with high accuracy (pp_data may slightly modify floating-point representation of number)."
    )
  }
})

check_inter_linpred_consistency <- function(example, model_data, logit_pref) {
  example$model_data <- model_data
  example$logit_pref <- logit_pref
  with(example, {
    fake_fit <- sample_prior_mean(blrmfit)
    drugs <- grep("^drug", names(model_data), value = TRUE)
    nr <- nrow(model_data)
    dcomp <- expand.grid(lapply(dref, c, 0))
    names(dcomp) <- drugs
    ## returns per term the log(1-p(DLT))
    drug_log_inverse_p_term <- function(term) {
      if (term == 0) {
        return(0)
      }
      return(log(inv_logit(-logit_pref)))
    }
    for (i in 1:nrow(dcomp)) {
      log_pNo <- sum(sapply(dcomp[i, ], drug_log_inverse_p_term))
      ref <- logit(1 - exp(log_pNo))
      nd <- model_data
      nd[, drugs] <- 0
      nd[, drugs] <- dcomp[i, ]
      pp <- posterior_linpred(fake_fit, newdata = nd, transform = FALSE)
      expect_equal(as.vector(pp[1, ]), rep(ref, times = nr))
    }
  })
}

test_that("posterior_linpred is consistent at reference dose (single-agent)", {
  single_agent <- load_single_agent()
  check_inter_linpred_consistency(single_agent, hist_SA, 0)
})
test_that("posterior_linpred is consistent at reference dose (combo2)", {
  combo2 <- load_combo2()
  check_inter_linpred_consistency(combo2, codata_combo2, logit(0.2))
})
test_that("posterior_linpred is consistent at reference dose (combo3)", {
  combo3 <- load_combo3()
  check_inter_linpred_consistency(combo3, hist_combo3, logit(1 / 3))
})

check_slope_linpred_consistency <- function(example, model_data, logit_pref) {
  example$model_data <- model_data
  example$logit_pref <- logit_pref
  with(example, {
    fake_fit <- sample_prior_mean(blrmfit)
    drugs <- grep("^drug", names(model_data), value = TRUE)
    nr <- nrow(model_data)
    dcomp <- expand.grid(lapply(exp(1) * dref, c, 0))
    names(dcomp) <- drugs
    ## returns per term the log(1-p(DLT))
    drug_log_inverse_p_term <- function(term) {
      if (term == 0) {
        return(0)
      }
      ## we assume that the prior mean for the slope is 1 here
      return(log(inv_logit(-logit_pref - 1)))
    }
    for (i in 1:nrow(dcomp)) {
      log_pNo <- sum(sapply(dcomp[i, ], drug_log_inverse_p_term))
      ref <- logit(1 - exp(log_pNo))
      nd <- model_data
      nd[, drugs] <- 0
      nd[, drugs] <- dcomp[i, ]
      pp <- posterior_linpred(fake_fit, newdata = nd, transform = FALSE)
      expect_equal(as.vector(pp[1, ]), rep(ref, times = nr))
    }
  })
}

test_that("posterior_linpred is consistent at exp(1) times reference dose (single-agent)", {
  single_agent <- load_single_agent()
  check_slope_linpred_consistency(single_agent, hist_SA, 0)
})
test_that("posterior_linpred is consistent at exp(1) times reference dose (combo2)", {
  combo2 <- load_combo2()
  check_slope_linpred_consistency(combo2, codata_combo2, logit(0.2))
})
test_that("posterior_linpred is consistent at exp(1) times reference dose (combo3)", {
  combo3 <- load_combo3()
  check_slope_linpred_consistency(combo3, hist_combo3, logit(1 / 3))
})

test_that("no unexpected error of posterior summary when num_groups = 1 or 2 and has_inter = TRUE (deprecated interface)", {
  skip_on_cran()

  withr::local_options(lifecycle_verbosity = "quiet")

  ## this test is to trigger a problem in posterior prior to
  ## 1.4.0, see https://github.com/stan-dev/posterior/issues/265
  ## for the bug report
  groups <- "one_group"
  combo_data <- tibble(
    group_id = factor(groups, groups),
    drug_A = 1,
    drug_B = 1,
    num_patients = 3,
    num_toxicities = 1
  )

  num_groups <- 1
  suppressWarnings(
    fit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug_A)) |
          1 + I(log(drug_B)) |
          0 + I(2 * drug_A * drug_B / (1 + drug_A * drug_B)) |
          group_id,
      data = combo_data,
      prior_EX_mu_mean_comp = matrix(
        c(
          log(1 / 4),
          0,
          log(1 / 4),
          0
        ),
        nrow = 2,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_mu_sd_comp = matrix(
        c(
          1,
          0.7,
          1,
          0.7
        ),
        nrow = 2,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_mu_mean_inter = 0,
      prior_EX_mu_sd_inter = 1,
      prior_EX_tau_mean_comp = matrix(
        c(0, 0),
        nrow = 2,
        ncol = 2
      ),
      prior_EX_tau_sd_comp = matrix(
        c(1, 1),
        nrow = 2,
        ncol = 2
      ),
      prior_EX_tau_mean_inter = matrix(0),
      prior_EX_tau_sd_inter = matrix(1),
      prior_is_EXNEX_comp = c(FALSE, FALSE),
      prior_is_EXNEX_inter = FALSE,
      prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = 2),
      prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = 1),
      prior_tau_dist = 0,
      prior_PD = FALSE
    )
  )

  expect_data_frame(summary(fit), nrows = 1)

  groups <- c(groups, "new_group")
  new_data <- combo_data
  levels(new_data$group_id) <- groups
  num_groups <- 2

  suppressWarnings(fit2 <- update(fit, data = new_data))
  expect_data_frame(summary(fit2), nrows = 1)
})


test_that("no unexpected error of posterior summary when num_groups = 1 or 2 and has_inter = TRUE", {
  skip_on_cran()

  ## this test is to trigger a problem in posterior prior to
  ## 1.4.0, see https://github.com/stan-dev/posterior/issues/265
  ## for the bug report
  groups <- "one_group"
  combo_data <- tibble(
    group_id = factor(groups, groups),
    drug_A = 1,
    drug_B = 1,
    num_patients = 3,
    num_toxicities = 1
  )

  num_groups <- 1
  num_comp <- 2
  num_inter <- 1

  suppressWarnings(
    fit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug_A)) |
          1 + I(log(drug_B)) |
          0 + I(2 * drug_A * drug_B / (1 + drug_A * drug_B)) |
          group_id,
      data = combo_data,
      prior_EX_mu_comp = list(
        mixmvnorm(c(1, log(1 / 4), 0, diag(c(1, 0.7)^2))),
        mixmvnorm(c(1, logit(1 / 4), 0, diag(c(1, 0.7)^2)))
      ),
      prior_EX_tau_comp = list(
        mixmvnorm(c(1, 0, 0, diag(c(1, 1)))),
        mixmvnorm(c(1, 0, 0, diag(c(1, 1))))
      ),
      prior_EX_mu_inter = mixmvnorm(c(1, 0, 1)),
      prior_EX_tau_inter = mixmvnorm(c(1, log(2) / 1.96, (log(2) / 1.96)^2)),
      prior_is_EXNEX_comp = rep(FALSE, num_comp),
      prior_is_EXNEX_inter = rep(FALSE, num_inter),
      prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = 2),
      prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = 1),
      prior_tau_dist = 0,
      init = 0,
      prior_PD = FALSE
    )
  )

  expect_data_frame(summary(fit), nrows = 1)

  groups <- c(groups, "new_group")
  new_data <- combo_data
  levels(new_data$group_id) <- groups
  num_groups <- 2

  suppressWarnings(fit2 <- update(fit, data = new_data))
  expect_data_frame(summary(fit2), nrows = 1)
})


## tests on prior consistency

test_that("specified prior of a dual combination with stratification matches sampled prior (deprecated interface)", {
  skip_on_cran()

  withr::local_options(lifecycle_verbosity = "quiet")
  withr::local_seed(67894356)

  hist_data <- tibble(
    group_id = as.factor(c(
      rep("trial_a", 2),
      rep("trial_b", 3),
      rep("trial_c", 1)
    )),
    stratum_id = as.factor(c(rep("reg1", 2), rep("reg2", 2), rep("reg2", 2))),
    drug1 = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 0),
    drug2 = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 10),
    num_toxicities = c(0, 1, 1, 0, 1, 0),
    num_patients = c(2, 6, 3, 4, 9, 29)
  )

  num_comp <- 2
  num_strata <- nlevels(hist_data$stratum_id)
  num_groups <- nlevels(hist_data$group_id)
  num_inter <- 1

  blrmfit_prior_PD <- blrm_exnex(
    cbind(num_toxicities, num_patients - num_toxicities) ~
      1 +
        I(log(drug1 / 100)) |
        1 + I(log(drug2 / 100)) |
        0 + I((drug1 / 100) * (drug2 / 100)) |
        stratum_id / group_id,
    data = hist_data,
    prior_is_EXNEX_comp = rep(FALSE, 2),
    prior_EX_prob_comp = matrix(
      c(1, 1, 1),
      nrow = num_groups,
      ncol = num_comp,
      byrow = FALSE
    ), # 0.5 would be ignored
    prior_is_EXNEX_inter = FALSE,
    prior_EX_prob_inter = matrix(
      c(1, 1, 1),
      nrow = num_groups,
      ncol = num_inter
    ),
    prior_EX_mu_mean_inter = rep(0, num_inter),
    prior_EX_mu_sd_inter = rep(1, num_inter),
    prior_EX_tau_mean_inter = matrix(
      log(2) / 1.96,
      nrow = num_strata,
      ncol = num_inter
    ),
    prior_EX_tau_sd_inter = matrix(
      log(2) / 1.96,
      nrow = num_strata,
      ncol = num_inter
    ),
    prior_EX_mu_mean_comp = matrix(
      c(logit(0.20), 0, logit(0.40), log(1.5)), # (E(mu_alpha), E(mu_beta))
      nrow = num_comp,
      ncol = 2,
      byrow = TRUE
    ),
    prior_EX_mu_sd_comp = matrix(
      c(1, 0.8, 2, 1) / 8, # (sd(mu_alpha), sd(mu_beta))
      nrow = num_comp,
      ncol = 2,
      byrow = TRUE
    ),
    prior_EX_tau_mean_comp = abind(
      matrix(log(c(0.25, 0.125) / 2), nrow = num_comp, ncol = 2, TRUE), # level 1 reg1
      matrix(log(2 * c(0.25, 0.125)), nrow = num_comp, ncol = 2, TRUE), # level 2 reg2
      along = 0
    ),
    prior_EX_tau_sd_comp = abind(
      matrix(
        c(log(2) / 1.96, log(1.25) / 1.96),
        nrow = num_comp,
        ncol = 2,
        TRUE
      ),
      matrix(
        c(log(3) / 1.96, log(2.5) / 1.96),
        nrow = num_comp,
        ncol = 2,
        TRUE
      ),
      along = 0
    ),
    prior_tau_dist = 1,
    prior_PD = TRUE,
    iter = 3000,
    warmup = 1000,
    chains = 1,
    cores = 1,
    init = 0,
    save_warmup = FALSE
  )

  blrmfit_data_is_like_prior <- update(
    blrmfit_prior_PD,
    data = mutate(hist_data, num_toxicities = 0, num_patients = 0),
    prior_PD = FALSE,
    iter = 3000,
    warmup = 1000,
    chains = 1,
    init = 0,
    save_warmup = FALSE
  )

  ## checks that the median of test is within 3 se of the referenc
  ## median and calculate a p-value from a ks test
  check_parameter <- function(test, ref) {
    param_name <- deparse(substitute(ref))
    test <- posterior::thin_draws(test, ceiling(ndraws(test) / 500))
    ref <- posterior::thin_draws(ref, ceiling(ndraws(ref) / 500))
    m_test <- unname(median(test))
    m_ref <- unname(median(ref))
    se_ref <- posterior::mcse_median(ref)
    se_test <- posterior::mcse_median(test)
    z <- (m_test - m_ref) / sqrt(se_test^2 + se_ref^2)
    expect_number(
      z,
      lower = qnorm(0.0005),
      upper = qnorm(0.9995),
      finite = TRUE,
      label = param_name,
      info = paste0("ref = ", ref, "; test = ", test, "; z = ", z)
    )
    ##suppressWarnings(ks.test(x=draws_of(test)[,1], y=draws_of(ref)[,1]))$p.value
    pnorm(-1 * abs(z))
  }

  check_sampled_prior <- function(blrmfit) {
    ref_data1 <- mutate(
      unique(hist_data[c("group_id", "stratum_id")]),
      num_toxicities = 0,
      num_patients = 0,
      drug1 = 100,
      drug2 = 0
    )
    ref_data2 <- mutate(
      unique(hist_data[c("group_id", "stratum_id")]),
      num_toxicities = 0,
      num_patients = 0,
      drug1 = 0,
      drug2 = 100
    )

    ## get group specific intercepts and slopes
    rv_prior_drug1_inter <- rvar(posterior_linpred(
      blrmfit,
      newdata = ref_data1
    ))
    rv_prior_drug1_log_slope <- log(
      rvar(posterior_linpred(
        blrmfit,
        newdata = mutate(ref_data1, drug1 = exp(1) * 100)
      )) -
        rv_prior_drug1_inter
    )

    rv_prior_drug2_inter <- rvar(posterior_linpred(
      blrmfit,
      newdata = ref_data2
    ))
    rv_prior_drug2_log_slope <- log(
      rvar(posterior_linpred(
        blrmfit,
        newdata = mutate(ref_data2, drug2 = exp(1) * 100)
      )) -
        rv_prior_drug2_inter
    )

    S <- ndraws(rv_prior_drug1_inter)

    prior_EX_mu_drug1_intercept <- rvar(rnorm(S, logit(0.2), 1 / 8))
    prior_EX_mu_drug1_log_slope <- rvar(rnorm(S, 0, 0.8 / 8))

    prior_EX_mu_drug2_intercept <- rvar(rnorm(S, logit(0.4), 2 / 8))
    prior_EX_mu_drug2_log_slope <- rvar(rnorm(S, log(1.5), 1 / 8))

    prior_EX_tau_intercept_strat1 <- rvar(exp(rnorm(
      S,
      log(0.250 / 2),
      log(2) / 1.96
    )))
    prior_EX_tau_log_slope_strat1 <- rvar(exp(rnorm(
      S,
      log(0.125 / 2),
      log(1.25) / 1.96
    )))

    prior_EX_tau_intercept_strat2 <- rvar(exp(rnorm(
      S,
      log(2 * 0.250),
      log(3) / 1.96
    )))
    prior_EX_tau_log_slope_strat2 <- rvar(exp(rnorm(
      S,
      log(2 * 0.125),
      log(2.5) / 1.96
    )))

    prior_drug1_intercept_strat1 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug1_intercept,
      sd = prior_EX_tau_intercept_strat1
    )
    prior_drug1_log_slope_strat1 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug1_log_slope,
      sd = prior_EX_tau_log_slope_strat1
    )

    prior_drug1_intercept_strat2 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug1_intercept,
      sd = prior_EX_tau_intercept_strat2
    )
    prior_drug1_log_slope_strat2 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug1_log_slope,
      sd = prior_EX_tau_log_slope_strat2
    )

    prior_drug2_intercept_strat1 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug2_intercept,
      sd = prior_EX_tau_intercept_strat1
    )
    prior_drug2_log_slope_strat1 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug2_log_slope,
      sd = prior_EX_tau_log_slope_strat1
    )

    prior_drug2_intercept_strat2 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug2_intercept,
      sd = prior_EX_tau_intercept_strat2
    )
    prior_drug2_log_slope_strat2 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug2_log_slope,
      sd = prior_EX_tau_log_slope_strat2
    )

    p_values <- c()
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_inter[1], prior_drug1_intercept_strat1)
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_log_slope[1], prior_drug1_log_slope_strat1)
    )

    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_inter[1], prior_drug2_intercept_strat1)
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_log_slope[1], prior_drug2_log_slope_strat1)
    )

    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_inter[2], prior_drug1_intercept_strat2)
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_log_slope[2], prior_drug1_log_slope_strat2)
    )

    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_inter[2], prior_drug2_intercept_strat2)
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_log_slope[2], prior_drug2_log_slope_strat2)
    )

    num_failed <- sum(p_values < 0.01)
    num_tests <- length(p_values)
    ## Aggregate calibration check on a stochastic MCMC-vs-prior
    ## comparison. Bound the count of nominal 1% deviations at the 0.1%
    ## quantile of its binomial null so this check tolerates the Monte
    ## Carlo noise of the sampled prior (notably the high-variance medians
    ## of the EXNEX mixture components), consistent with the per-parameter
    ## hard bound on |z| enforced above.
    expect_integer(
      num_failed,
      lower = 0,
      upper = qbinom(0.999, num_tests, 0.01),
      any.missing = FALSE
    )
  }

  check_sampled_prior(blrmfit_prior_PD)
  check_sampled_prior(blrmfit_data_is_like_prior)
})


test_that("specified prior of a dual combination with stratification matches sampled prior", {
  skip_on_cran()

  withr::local_seed(67894356)

  hist_data <- tibble(
    group_id = as.factor(c(
      rep("trial_a", 2),
      rep("trial_b", 3),
      rep("trial_c", 1)
    )),
    stratum_id = as.factor(c(rep("reg1", 2), rep("reg2", 2), rep("reg2", 2))),
    drug1 = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 0),
    drug2 = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 10),
    num_toxicities = c(0, 1, 1, 0, 1, 0),
    num_patients = c(2, 6, 3, 4, 9, 29)
  )

  num_comp <- 2
  num_strata <- nlevels(hist_data$stratum_id)
  num_groups <- nlevels(hist_data$group_id)
  num_inter <- 1

  blrmfit_prior_PD <- blrm_exnex(
    cbind(num_toxicities, num_patients - num_toxicities) ~
      1 +
        I(log(drug1 / 100)) |
        1 + I(log(drug2 / 100)) |
        0 + I(drug1 / 100 * drug2 / 100) |
        stratum_id / group_id,
    data = hist_data,
    prior_EX_mu_comp = list(
      mixmvnorm(c(1, logit(0.2), 0, diag(c(1 / 8, 0.8 / 8)^2))),
      mixmvnorm(c(1, logit(0.4), log(1.5), diag(c(2 / 8, 1 / 8)^2)))
    ),
    prior_EX_tau_comp = list(
      list(
        mixmvnorm(c(
          1,
          log(0.250 / 2),
          log(0.125 / 2),
          diag(c(log(2) / 1.96, log(1.25) / 1.96)^2)
        )),
        mixmvnorm(c(
          1,
          log(0.250 / 2),
          log(0.125 / 2),
          diag(c(log(2) / 1.96, log(1.25) / 1.96)^2)
        ))
      ),
      list(
        mixmvnorm(c(
          1,
          log(2 * 0.250),
          log(2 * 0.125),
          diag(c(log(3) / 1.96, log(2.5) / 1.96)^2)
        )),
        mixmvnorm(c(
          1,
          log(2 * 0.250),
          log(2 * 0.125),
          diag(c(log(3) / 1.96, log(2.5) / 1.96)^2)
        ))
      )
    ),
    prior_EX_mu_inter = mixmvnorm(c(1, 0, 1)),
    prior_EX_tau_inter = list(
      mixmvnorm(c(1, log(2) / 1.96, (log(2) / 1.96)^2)),
      mixmvnorm(c(1, log(2) / 1.96, (log(2) / 1.96)^2))
    ),
    prior_is_EXNEX_comp = rep(FALSE, num_comp),
    prior_is_EXNEX_inter = rep(FALSE, num_inter),
    prior_EX_prob_comp = cbind(c(1, 1, 1), c(1, 1, 1)),
    prior_EX_prob_inter = matrix(
      c(1, 1, 1),
      nrow = num_groups,
      ncol = num_inter
    ),
    prior_tau_dist = 1,
    prior_PD = TRUE,
    iter = 3000,
    warmup = 1000,
    chains = 1,
    init = 0,
    save_warmup = FALSE
  )

  blrmfit_data_is_like_prior <- update(
    blrmfit_prior_PD,
    data = mutate(hist_data, num_toxicities = 0, num_patients = 0),
    prior_PD = FALSE,
    iter = 3000,
    warmup = 1000,
    chains = 1,
    init = 0,
    save_warmup = FALSE
  )

  ## checks that the median of test is within 3 se of the referenc
  ## median and calculate a p-value from a ks test
  check_parameter <- function(test, ref) {
    param_name <- deparse(substitute(ref))
    test <- posterior::thin_draws(test, ceiling(ndraws(test) / 500))
    ref <- posterior::thin_draws(ref, ceiling(ndraws(ref) / 500))
    m_test <- unname(median(test))
    m_ref <- unname(median(ref))
    se_ref <- posterior::mcse_median(ref)
    se_test <- posterior::mcse_median(test)
    z <- (m_test - m_ref) / sqrt(se_test^2 + se_ref^2)
    expect_number(
      z,
      lower = qnorm(0.0005),
      upper = qnorm(0.9995),
      finite = TRUE,
      label = param_name,
      info = paste0("ref = ", ref, "; test = ", test, "; z = ", z)
    )
    ##suppressWarnings(ks.test(x=draws_of(test)[,1], y=draws_of(ref)[,1]))$p.value
    pnorm(-1 * abs(z))
  }

  check_sampled_prior <- function(blrmfit) {
    ref_data1 <- mutate(
      unique(hist_data[c("group_id", "stratum_id")]),
      num_toxicities = 0,
      num_patients = 0,
      drug1 = 100,
      drug2 = 0
    )
    ref_data2 <- mutate(
      unique(hist_data[c("group_id", "stratum_id")]),
      num_toxicities = 0,
      num_patients = 0,
      drug1 = 0,
      drug2 = 100
    )

    ## get group specific intercepts and slopes
    rv_prior_drug1_inter <- rvar(posterior_linpred(
      blrmfit,
      newdata = ref_data1
    ))
    rv_prior_drug1_log_slope <- log(
      rvar(posterior_linpred(
        blrmfit,
        newdata = mutate(ref_data1, drug1 = exp(1) * 100)
      )) -
        rv_prior_drug1_inter
    )

    rv_prior_drug2_inter <- rvar(posterior_linpred(
      blrmfit,
      newdata = ref_data2
    ))
    rv_prior_drug2_log_slope <- log(
      rvar(posterior_linpred(
        blrmfit,
        newdata = mutate(ref_data2, drug2 = exp(1) * 100)
      )) -
        rv_prior_drug2_inter
    )

    S <- ndraws(rv_prior_drug1_inter)

    prior_EX_mu_drug1_intercept <- rvar(rnorm(S, logit(0.2), 1 / 8))
    prior_EX_mu_drug1_log_slope <- rvar(rnorm(S, 0, 0.8 / 8))

    prior_EX_mu_drug2_intercept <- rvar(rnorm(S, logit(0.4), 2 / 8))
    prior_EX_mu_drug2_log_slope <- rvar(rnorm(S, log(1.5), 1 / 8))

    prior_EX_tau_intercept_strat1 <- rvar(exp(rnorm(
      S,
      log(0.250 / 2),
      log(2) / 1.96
    )))
    prior_EX_tau_log_slope_strat1 <- rvar(exp(rnorm(
      S,
      log(0.125 / 2),
      log(1.25) / 1.96
    )))

    prior_EX_tau_intercept_strat2 <- rvar(exp(rnorm(
      S,
      log(2 * 0.250),
      log(3) / 1.96
    )))
    prior_EX_tau_log_slope_strat2 <- rvar(exp(rnorm(
      S,
      log(2 * 0.125),
      log(2.5) / 1.96
    )))

    prior_drug1_intercept_strat1 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug1_intercept,
      sd = prior_EX_tau_intercept_strat1
    )
    prior_drug1_log_slope_strat1 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug1_log_slope,
      sd = prior_EX_tau_log_slope_strat1
    )

    prior_drug1_intercept_strat2 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug1_intercept,
      sd = prior_EX_tau_intercept_strat2
    )
    prior_drug1_log_slope_strat2 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug1_log_slope,
      sd = prior_EX_tau_log_slope_strat2
    )

    prior_drug2_intercept_strat1 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug2_intercept,
      sd = prior_EX_tau_intercept_strat1
    )
    prior_drug2_log_slope_strat1 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug2_log_slope,
      sd = prior_EX_tau_log_slope_strat1
    )

    prior_drug2_intercept_strat2 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug2_intercept,
      sd = prior_EX_tau_intercept_strat2
    )
    prior_drug2_log_slope_strat2 <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug2_log_slope,
      sd = prior_EX_tau_log_slope_strat2
    )

    p_values <- c()
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_inter[1], prior_drug1_intercept_strat1)
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_log_slope[1], prior_drug1_log_slope_strat1)
    )

    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_inter[1], prior_drug2_intercept_strat1)
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_log_slope[1], prior_drug2_log_slope_strat1)
    )

    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_inter[2], prior_drug1_intercept_strat2)
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_log_slope[2], prior_drug1_log_slope_strat2)
    )

    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_inter[2], prior_drug2_intercept_strat2)
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_log_slope[2], prior_drug2_log_slope_strat2)
    )

    num_failed <- sum(p_values < 0.01)
    num_tests <- length(p_values)
    ## Aggregate calibration check on a stochastic MCMC-vs-prior
    ## comparison. Bound the count of nominal 1% deviations at the 0.1%
    ## quantile of its binomial null so this check tolerates the Monte
    ## Carlo noise of the sampled prior (notably the high-variance medians
    ## of the EXNEX mixture components), consistent with the per-parameter
    ## hard bound on |z| enforced above.
    expect_integer(
      num_failed,
      lower = 0,
      upper = qbinom(0.999, num_tests, 0.01),
      any.missing = FALSE
    )
  }

  check_sampled_prior(blrmfit_prior_PD)
  check_sampled_prior(blrmfit_data_is_like_prior)
})

test_that("specified prior of a dual combination with EXNEX matches sampled prior (deprecated interface)", {
  skip_on_cran()
  withr::local_seed(67894356)
  withr::local_options(lifecycle_verbosity = "quiet")

  hist_data <- tibble(
    group_id = as.factor(c(
      rep("trial_a", 2),
      rep("trial_b", 3),
      rep("trial_c", 1)
    )),
    stratum_id = as.factor(c(rep("reg1", 2), rep("reg1", 2), rep("reg1", 2))),
    drug1 = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 0),
    drug2 = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 10),
    num_toxicities = c(0, 1, 1, 0, 1, 0),
    num_patients = c(2, 6, 3, 4, 9, 29)
  )

  num_comp <- 2
  num_strata <- nlevels(hist_data$stratum_id)
  num_groups <- nlevels(hist_data$group_id)
  num_inter <- 1

  blrmfit_prior_PD <- blrm_exnex(
    cbind(num_toxicities, num_patients - num_toxicities) ~
      1 +
        I(log(drug1 / 100)) |
        1 + I(log(drug2 / 100)) |
        0 + I((drug1 / 100) * (drug2 / 100)) |
        stratum_id / group_id,
    data = hist_data,
    prior_is_EXNEX_comp = rep(TRUE, 2),
    prior_EX_prob_comp = matrix(
      c(1, 0.5, 0.1, 0.2, 0.5, 1),
      nrow = num_groups,
      ncol = num_comp,
      byrow = FALSE
    ),
    prior_is_EXNEX_inter = FALSE,
    prior_EX_prob_inter = matrix(
      c(1, 1, 1),
      nrow = num_groups,
      ncol = num_inter
    ),
    prior_EX_mu_mean_inter = rep(0, num_inter),
    prior_EX_mu_sd_inter = rep(1, num_inter),
    prior_EX_tau_mean_inter = matrix(
      log(2) / 1.96,
      nrow = num_strata,
      ncol = num_inter
    ),
    prior_EX_tau_sd_inter = matrix(
      log(2) / 1.96,
      nrow = num_strata,
      ncol = num_inter
    ),
    prior_EX_mu_mean_comp = matrix(
      c(logit(0.20), 0, logit(0.40), log(1.5)), # (E(mu_alpha), E(mu_beta))
      nrow = num_comp,
      ncol = 2,
      byrow = TRUE
    ),
    prior_EX_mu_sd_comp = matrix(
      c(1, 0.5, 2, 1) / 8, # (sd(mu_alpha), sd(mu_beta))
      nrow = num_comp,
      ncol = 2,
      byrow = TRUE
    ),
    prior_EX_tau_mean_comp = abind(
      matrix(log(c(0.25, 0.125) / 2), nrow = num_comp, ncol = 2, TRUE), # level 1 reg1
      along = 0
    ),
    prior_EX_tau_sd_comp = abind(
      matrix(
        c(log(2) / 1.96, log(1.25) / 1.96),
        nrow = num_comp,
        ncol = 2,
        TRUE
      ),
      along = 0
    ),
    prior_NEX_mu_mean_comp = matrix(
      c(logit(c(2 / 3, 3 / 4)), 0.2, 0.4),
      nrow = num_comp,
      ncol = 2,
      byrow = FALSE
    ),
    prior_NEX_mu_sd_comp = matrix(
      c(0.5, 0.25, 0.2, 0.1),
      nrow = num_comp,
      ncol = 2,
      byrow = FALSE
    ),
    prior_tau_dist = 1,
    prior_PD = TRUE,
    iter = 4000,
    warmup = 2000,
    chains = 2,
    cores = 1,
    init = 0,
    save_warmup = FALSE
  )

  blrmfit_data_is_prior_like <- update(
    blrmfit_prior_PD,
    data = mutate(hist_data, num_toxicities = 0, num_patients = 0),
    prior_PD = FALSE,
    iter = 4000,
    warmup = 2000,
    chains = 2,
    init = 0,
    save_warmup = FALSE
  )

  rv_sample <- function(s1, s2, w1)
    rvar_ifelse(
      rvar(sample.int(2, ndraws(s1), replace = TRUE, prob = c(w1, 1 - w1))) ==
        1,
      s1,
      s2
    )

  ref_data1 <- mutate(
    unique(hist_data[c("group_id", "stratum_id")]),
    num_toxicities = 0,
    num_patients = 0,
    drug1 = 100,
    drug2 = 0
  )
  ref_data2 <- mutate(
    unique(hist_data[c("group_id", "stratum_id")]),
    num_toxicities = 0,
    num_patients = 0,
    drug1 = 0,
    drug2 = 100
  )

  ## checks that the median of test is within 3 se of the reference
  check_parameter <- function(test, ref) {
    param_name <- deparse(substitute(ref))
    test <- posterior::thin_draws(test, ceiling(ndraws(test) / 1000))
    ref <- posterior::thin_draws(ref, ceiling(ndraws(ref) / 1000))
    m_test <- unname(median(test))
    m_ref <- unname(median(ref))
    se_ref <- posterior::mcse_median(ref)
    se_test <- posterior::mcse_median(test)
    z <- (m_test - m_ref) / sqrt(se_test^2 + se_ref^2)
    expect_number(
      z,
      lower = -6,
      upper = 6,
      finite = TRUE,
      label = param_name,
      info = paste0("ref = ", ref, "; test = ", test, "; z = ", z)
    )
    ##suppressWarnings(ks.test(x=draws_of(test)[,1], y=draws_of(ref)[,1]))$p.value
    pnorm(-1 * abs(z))
  }

  check_sampled_prior <- function(blrmfit) {
    ## get group specific intercepts and slopes
    rv_prior_drug1_inter <- rvar(posterior_linpred(
      blrmfit,
      newdata = ref_data1
    ))
    rv_prior_drug1_log_slope <- log(
      rvar(posterior_linpred(
        blrmfit,
        newdata = mutate(ref_data1, drug1 = exp(1) * 100)
      )) -
        rv_prior_drug1_inter
    )

    rv_prior_drug2_inter <- rvar(posterior_linpred(
      blrmfit,
      newdata = ref_data2
    ))
    rv_prior_drug2_log_slope <- log(
      rvar(posterior_linpred(
        blrmfit,
        newdata = mutate(ref_data2, drug2 = exp(1) * 100)
      )) -
        rv_prior_drug2_inter
    )

    S <- ndraws(rv_prior_drug1_inter)

    prior_EX_mu_drug1_intercept <- rvar(rnorm(S, logit(0.2), 1 / 8))
    prior_EX_mu_drug1_log_slope <- rvar(rnorm(S, 0, 0.5 / 8))

    prior_EX_mu_drug2_intercept <- rvar(rnorm(S, logit(0.4), 2 / 8))
    prior_EX_mu_drug2_log_slope <- rvar(rnorm(S, log(1.5), 1 / 8))

    prior_EX_tau_intercept <- rvar(exp(rnorm(S, log(0.250 / 2), log(2) / 1.96)))
    prior_EX_tau_log_slope <- rvar(exp(rnorm(
      S,
      log(0.125 / 2),
      log(1.25) / 1.96
    )))

    prior_NEX_drug1_intercept <- rvar(rnorm(S, logit(2 / 3), 0.5))
    prior_NEX_drug1_log_slope <- rvar(rnorm(S, 0.2, 0.25))

    prior_NEX_drug2_intercept <- rvar(rnorm(S, logit(3 / 4), 0.2))
    prior_NEX_drug2_log_slope <- rvar(rnorm(S, 0.4, 0.1))

    prior_EX_drug1_intercept <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug1_intercept,
      sd = prior_EX_tau_intercept
    )
    prior_EX_drug1_log_slope <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug1_log_slope,
      sd = prior_EX_tau_log_slope
    )

    prior_EX_drug2_intercept <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug2_intercept,
      sd = prior_EX_tau_intercept
    )
    prior_EX_drug2_log_slope <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug2_log_slope,
      sd = prior_EX_tau_log_slope
    )

    p_values <- c()

    prior_EXNEX_drug1_intercept <- rv_sample(
      prior_EX_drug1_intercept,
      prior_NEX_drug1_intercept,
      1
    )
    ## abs(rv_prior_drug1_inter[1] - prior_EXNEX_drug1_intercept)
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_inter[1], prior_EXNEX_drug1_intercept)
    )

    prior_EXNEX_drug1_log_slope <- rv_sample(
      prior_EX_drug1_log_slope,
      prior_NEX_drug1_log_slope,
      1
    )
    ## abs(rv_prior_drug1_log_slope[1] - prior_EXNEX_drug1_log_slope)
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_log_slope[1], prior_EXNEX_drug1_log_slope)
    )

    prior_EXNEX_drug2_intercept <- rv_sample(
      prior_EX_drug2_intercept,
      prior_NEX_drug2_intercept,
      0.2
    )
    ##abs(rv_prior_drug2_inter[1] - prior_EXNEX_drug2_intercept)
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_inter[1], prior_EXNEX_drug2_intercept)
    )

    prior_EXNEX_drug2_log_slope <- rv_sample(
      prior_EX_drug2_log_slope,
      prior_NEX_drug2_log_slope,
      0.2
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_log_slope[1], prior_EXNEX_drug2_log_slope)
    )

    prior_EXNEX_drug1_intercept <- rv_sample(
      prior_EX_drug1_intercept,
      prior_NEX_drug1_intercept,
      0.5
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_inter[2], prior_EXNEX_drug1_intercept)
    )

    prior_EXNEX_drug1_log_slope <- rv_sample(
      prior_EX_drug1_log_slope,
      prior_NEX_drug1_log_slope,
      0.5
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_log_slope[2], prior_EXNEX_drug1_log_slope)
    )

    prior_EXNEX_drug2_intercept <- rv_sample(
      prior_EX_drug2_intercept,
      prior_NEX_drug2_intercept,
      0.5
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_inter[2], prior_EXNEX_drug2_intercept)
    )

    prior_EXNEX_drug2_log_slope <- rv_sample(
      prior_EX_drug2_log_slope,
      prior_NEX_drug2_log_slope,
      0.5
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_log_slope[2], prior_EXNEX_drug2_log_slope)
    )

    prior_EXNEX_drug1_intercept <- rv_sample(
      prior_EX_drug1_intercept,
      prior_NEX_drug1_intercept,
      0.1
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_inter[3], prior_EXNEX_drug1_intercept)
    )

    prior_EXNEX_drug1_log_slope <- rv_sample(
      prior_EX_drug1_log_slope,
      prior_NEX_drug1_log_slope,
      0.1
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_log_slope[3], prior_EXNEX_drug1_log_slope)
    )

    prior_EXNEX_drug2_intercept <- rv_sample(
      prior_EX_drug2_intercept,
      prior_NEX_drug2_intercept,
      1.0
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_inter[3], prior_EXNEX_drug2_intercept)
    )

    prior_EXNEX_drug2_log_slope <- rv_sample(
      prior_EX_drug2_log_slope,
      prior_NEX_drug2_log_slope,
      1.0
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_log_slope[3], prior_EXNEX_drug2_log_slope)
    )

    num_failed <- sum(p_values < 0.01)
    num_tests <- length(p_values)
    ## Aggregate calibration check on a stochastic MCMC-vs-prior
    ## comparison. Bound the count of nominal 1% deviations at the 0.1%
    ## quantile of its binomial null so this check tolerates the Monte
    ## Carlo noise of the sampled prior (notably the high-variance medians
    ## of the EXNEX mixture components), consistent with the per-parameter
    ## hard bound on |z| enforced above.
    expect_integer(
      num_failed,
      lower = 0,
      upper = qbinom(0.999, num_tests, 0.01),
      any.missing = FALSE
    )
  }

  check_sampled_prior(blrmfit_prior_PD)
  check_sampled_prior(blrmfit_data_is_prior_like)
})


test_that("specified prior of a dual combination with EXNEX matches sampled prior", {
  skip_on_cran()
  withr::local_seed(67894356)

  hist_data <- tibble(
    group_id = as.factor(c(
      rep("trial_a", 2),
      rep("trial_b", 3),
      rep("trial_c", 1)
    )),
    stratum_id = as.factor(c(rep("reg1", 2), rep("reg1", 2), rep("reg1", 2))),
    drug1 = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 0),
    drug2 = c(20 * 5, 30 * 5, 20 * 14, 30 * 14, 45 * 7, 10),
    num_toxicities = c(0, 1, 1, 0, 1, 0),
    num_patients = c(2, 6, 3, 4, 9, 29)
  )

  num_comp <- 2
  num_strata <- nlevels(hist_data$stratum_id)
  num_groups <- nlevels(hist_data$group_id)
  num_inter <- 1

  blrmfit_prior_PD <- blrm_exnex(
    cbind(num_toxicities, num_patients - num_toxicities) ~
      1 +
        I(log(drug1 / 100)) |
        1 + I(log(drug2 / 100)) |
        0 + I(drug1 / 100 * drug2 / 100) |
        stratum_id / group_id,
    data = hist_data,
    prior_EX_mu_comp = list(
      mixmvnorm(c(1, logit(0.2), 0, diag(c(1 / 8, 0.5 / 8)^2))),
      mixmvnorm(c(1, logit(0.4), log(1.5), diag(c(2 / 8, 1 / 8)^2)))
    ),
    prior_EX_tau_comp = list(list(
      mixmvnorm(c(
        1,
        log(0.250 / 2),
        log(0.125 / 2),
        diag(c(log(2) / 1.96, log(1.25) / 1.96)^2)
      )),
      mixmvnorm(c(
        1,
        log(0.250 / 2),
        log(0.125 / 2),
        diag(c(log(2) / 1.96, log(1.25) / 1.96)^2)
      ))
    )),
    prior_EX_mu_inter = mixmvnorm(c(1, 0, 1)),
    prior_EX_tau_inter = list(mixmvnorm(c(
      1,
      log(2) / 1.96,
      (log(2) / 1.96)^2
    ))),
    prior_is_EXNEX_comp = rep(TRUE, num_comp),
    prior_EX_prob_comp = matrix(
      c(1, 0.5, 0.1, 0.2, 0.5, 1),
      nrow = num_groups,
      ncol = num_comp,
      byrow = FALSE
    ),
    prior_is_EXNEX_inter = rep(FALSE, num_inter),
    prior_EX_prob_inter = matrix(
      c(1, 1, 1),
      nrow = num_groups,
      ncol = num_inter
    ),
    prior_NEX_mu_comp = list(
      mixmvnorm(c(1, logit(2 / 3), 0.2, diag(c(0.5, 0.2)^2))),
      mixmvnorm(c(1, logit(3 / 4), 0.4, diag(c(0.25, 0.1)^2)))
    ),
    prior_tau_dist = 1,
    prior_PD = TRUE,
    iter = 4000,
    warmup = 2000,
    chains = 2,
    init = 0,
    save_warmup = FALSE
  )

  blrmfit_data_is_like_prior <- update(
    blrmfit_prior_PD,
    data = mutate(hist_data, num_toxicities = 0, num_patients = 0),
    prior_PD = FALSE,
    iter = 4000,
    warmup = 2000,
    chains = 2,
    init = 0,
    save_warmup = FALSE
  )

  rv_sample <- function(s1, s2, w1)
    rvar_ifelse(
      rvar(sample.int(2, ndraws(s1), replace = TRUE, prob = c(w1, 1 - w1))) ==
        1,
      s1,
      s2
    )

  ref_data1 <- mutate(
    unique(hist_data[c("group_id", "stratum_id")]),
    num_toxicities = 0,
    num_patients = 0,
    drug1 = 100,
    drug2 = 0
  )
  ref_data2 <- mutate(
    unique(hist_data[c("group_id", "stratum_id")]),
    num_toxicities = 0,
    num_patients = 0,
    drug1 = 0,
    drug2 = 100
  )

  ## checks that the median of test is within 3 se of the referenc
  ## median and calculate a p-value from a ks test
  check_parameter <- function(test, ref) {
    param_name <- deparse(substitute(ref))
    test <- posterior::thin_draws(test, ceiling(ndraws(test) / 1000))
    ref <- posterior::thin_draws(ref, ceiling(ndraws(ref) / 1000))
    m_test <- unname(median(test))
    m_ref <- unname(median(ref))
    se_ref <- posterior::mcse_median(ref)
    se_test <- posterior::mcse_median(test)
    z <- (m_test - m_ref) / sqrt(se_test^2 + se_ref^2)
    expect_number(
      z,
      lower = -6,
      upper = 6,
      finite = TRUE,
      label = param_name,
      info = paste0("ref = ", ref, "; test = ", test, "; z = ", z)
    )
    ##suppressWarnings(ks.test(x=draws_of(test)[,1], y=draws_of(ref)[,1]))$p.value
    pnorm(-1 * abs(z))
  }

  check_sampled_prior <- function(blrmfit) {
    ## get group specific intercepts and slopes
    rv_prior_drug1_inter <- rvar(posterior_linpred(
      blrmfit,
      newdata = ref_data1
    ))
    rv_prior_drug1_log_slope <- log(
      rvar(posterior_linpred(
        blrmfit,
        newdata = mutate(ref_data1, drug1 = exp(1) * 100)
      )) -
        rv_prior_drug1_inter
    )

    rv_prior_drug2_inter <- rvar(posterior_linpred(
      blrmfit,
      newdata = ref_data2
    ))
    rv_prior_drug2_log_slope <- log(
      rvar(posterior_linpred(
        blrmfit,
        newdata = mutate(ref_data2, drug2 = exp(1) * 100)
      )) -
        rv_prior_drug2_inter
    )

    S <- ndraws(rv_prior_drug1_inter)

    prior_EX_mu_drug1_intercept <- rvar(rnorm(S, logit(0.2), 1 / 8))
    prior_EX_mu_drug1_log_slope <- rvar(rnorm(S, 0, 0.5 / 8))

    prior_EX_mu_drug2_intercept <- rvar(rnorm(S, logit(0.4), 2 / 8))
    prior_EX_mu_drug2_log_slope <- rvar(rnorm(S, log(1.5), 1 / 8))

    prior_EX_tau_intercept <- rvar(exp(rnorm(S, log(0.250 / 2), log(2) / 1.96)))
    prior_EX_tau_log_slope <- rvar(exp(rnorm(
      S,
      log(0.125 / 2),
      log(1.25) / 1.96
    )))

    prior_NEX_drug1_intercept <- rvar(rnorm(S, logit(2 / 3), 0.5))
    prior_NEX_drug1_log_slope <- rvar(rnorm(S, 0.2, 0.25))

    prior_NEX_drug2_intercept <- rvar(rnorm(S, logit(3 / 4), 0.2))
    prior_NEX_drug2_log_slope <- rvar(rnorm(S, 0.4, 0.1))

    prior_EX_drug1_intercept <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug1_intercept,
      sd = prior_EX_tau_intercept
    )
    prior_EX_drug1_log_slope <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug1_log_slope,
      sd = prior_EX_tau_log_slope
    )

    prior_EX_drug2_intercept <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug2_intercept,
      sd = prior_EX_tau_intercept
    )
    prior_EX_drug2_log_slope <- rvar_rng(
      rnorm,
      1,
      mean = prior_EX_mu_drug2_log_slope,
      sd = prior_EX_tau_log_slope
    )

    p_values <- c()

    prior_EXNEX_drug1_intercept <- rv_sample(
      prior_EX_drug1_intercept,
      prior_NEX_drug1_intercept,
      1
    )
    ## abs(rv_prior_drug1_inter[1] - prior_EXNEX_drug1_intercept)
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_inter[1], prior_EXNEX_drug1_intercept)
    )

    prior_EXNEX_drug1_log_slope <- rv_sample(
      prior_EX_drug1_log_slope,
      prior_NEX_drug1_log_slope,
      1
    )
    ## abs(rv_prior_drug1_log_slope[1] - prior_EXNEX_drug1_log_slope)
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_log_slope[1], prior_EXNEX_drug1_log_slope)
    )

    prior_EXNEX_drug2_intercept <- rv_sample(
      prior_EX_drug2_intercept,
      prior_NEX_drug2_intercept,
      0.2
    )
    ##abs(rv_prior_drug2_inter[1] - prior_EXNEX_drug2_intercept)
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_inter[1], prior_EXNEX_drug2_intercept)
    )

    prior_EXNEX_drug2_log_slope <- rv_sample(
      prior_EX_drug2_log_slope,
      prior_NEX_drug2_log_slope,
      0.2
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_log_slope[1], prior_EXNEX_drug2_log_slope)
    )

    prior_EXNEX_drug1_intercept <- rv_sample(
      prior_EX_drug1_intercept,
      prior_NEX_drug1_intercept,
      0.5
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_inter[2], prior_EXNEX_drug1_intercept)
    )

    prior_EXNEX_drug1_log_slope <- rv_sample(
      prior_EX_drug1_log_slope,
      prior_NEX_drug1_log_slope,
      0.5
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_log_slope[2], prior_EXNEX_drug1_log_slope)
    )

    prior_EXNEX_drug2_intercept <- rv_sample(
      prior_EX_drug2_intercept,
      prior_NEX_drug2_intercept,
      0.5
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_inter[2], prior_EXNEX_drug2_intercept)
    )

    prior_EXNEX_drug2_log_slope <- rv_sample(
      prior_EX_drug2_log_slope,
      prior_NEX_drug2_log_slope,
      0.5
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_log_slope[2], prior_EXNEX_drug2_log_slope)
    )

    prior_EXNEX_drug1_intercept <- rv_sample(
      prior_EX_drug1_intercept,
      prior_NEX_drug1_intercept,
      0.1
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_inter[3], prior_EXNEX_drug1_intercept)
    )

    prior_EXNEX_drug1_log_slope <- rv_sample(
      prior_EX_drug1_log_slope,
      prior_NEX_drug1_log_slope,
      0.1
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug1_log_slope[3], prior_EXNEX_drug1_log_slope)
    )

    prior_EXNEX_drug2_intercept <- rv_sample(
      prior_EX_drug2_intercept,
      prior_NEX_drug2_intercept,
      1.0
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_inter[3], prior_EXNEX_drug2_intercept)
    )

    prior_EXNEX_drug2_log_slope <- rv_sample(
      prior_EX_drug2_log_slope,
      prior_NEX_drug2_log_slope,
      1.0
    )
    p_values <- c(
      p_values,
      check_parameter(rv_prior_drug2_log_slope[3], prior_EXNEX_drug2_log_slope)
    )

    num_failed <- sum(p_values < 0.01)
    num_tests <- length(p_values)
    ## Aggregate calibration check on a stochastic MCMC-vs-prior
    ## comparison. Bound the count of nominal 1% deviations at the 0.1%
    ## quantile of its binomial null so this check tolerates the Monte
    ## Carlo noise of the sampled prior (notably the high-variance medians
    ## of the EXNEX mixture components), consistent with the per-parameter
    ## hard bound on |z| enforced above.
    expect_integer(
      num_failed,
      lower = 0,
      upper = qbinom(0.999, num_tests, 0.01),
      any.missing = FALSE
    )
  }

  check_sampled_prior(blrmfit_prior_PD)
  check_sampled_prior(blrmfit_data_is_like_prior)
})

## test that we can call the old function arguments, but get a warning

test_that("blrm_exnex supports deprecated non-mixture arguments with a warning", {
  skip_on_cran()

  withr::local_seed(67894356)
  withr::local_options(lifecycle_verbosity = "quiet")

  dref <- c(6, 960)

  num_comp <- 2 # two investigational drugs
  num_inter <- 1 # one drug-drug interaction needs to be modeled
  num_groups <- nlevels(codata_combo2$group_id) # no stratification needed
  num_strata <- 1 # no stratification needed

  lifecycle::expect_deprecated(blrm_exnex(
    cbind(num_toxicities, num_patients - num_toxicities) ~
      1 +
        I(log(drug_A / dref[1])) |
        1 + I(log(drug_B / dref[2])) |
        0 + I(drug_A / dref[1] * drug_B / dref[2]) |
        group_id,
    data = codata_combo2,
    prior_EX_mu_mean_comp = matrix(
      c(
        logit(0.2),
        0, # hyper-mean of (intercept, log-slope) for drug A
        logit(0.2),
        0
      ), # hyper-mean of (intercept, log-slope) for drug B
      nrow = num_comp,
      ncol = 2,
      byrow = TRUE
    ),
    prior_EX_mu_sd_comp = matrix(
      c(
        2.0,
        1, # hyper-sd of mean mu for (intercept, log-slope) for drug A
        2.0,
        1
      ), # hyper-sd of mean mu for (intercept, log-slope) for drug B
      nrow = num_comp,
      ncol = 2,
      byrow = TRUE
    ),
    prior_EX_tau_mean_comp = matrix(
      0,
      nrow = num_comp,
      ncol = 2,
      byrow = TRUE
    ),
    prior_EX_tau_sd_comp = matrix(
      c(log(4) / 1.96, log(4) / 1.96, log(4) / 1.96, log(4) / 1.96),
      nrow = num_comp,
      ncol = 2,
      byrow = TRUE
    ),
    prior_EX_mu_mean_inter = 0,
    prior_EX_mu_sd_inter = 1.121,
    prior_EX_tau_mean_inter = matrix(0, nrow = num_strata, ncol = num_inter),
    prior_EX_tau_sd_inter = matrix(
      log(4) / 1.96,
      nrow = num_strata,
      ncol = num_inter
    ),
    prior_is_EXNEX_comp = rep(FALSE, num_comp),
    prior_is_EXNEX_inter = rep(FALSE, num_inter),
    prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = num_comp),
    prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = num_inter),
    ## setup sampling to be fast, but do not trigger any rstan warning
    ## as the rstan warning would otherwise
    prior_tau_dist = 0,
    prior_PD = TRUE,
    init = 0,
    chains = 1,
    iter = 1000,
    warmup = 500,
    control = list(adapt_delta = 0.85)
  ))
})

test_that("blrm_exnex does allow either non-mixture or mixture arguments, but no mix", {
  skip_on_cran()

  withr::local_seed(67894356)
  withr::local_options(lifecycle_verbosity = "quiet")

  dref <- c(6, 960)

  num_comp <- 2 # two investigational drugs
  num_inter <- 1 # one drug-drug interaction needs to be modeled
  num_groups <- nlevels(codata_combo2$group_id) # no stratification needed
  num_strata <- 1 # no stratification needed

  expect_error(
    suppressWarnings(blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug_A / dref[1])) |
          1 + I(log(drug_B / dref[2])) |
          0 + I(drug_A / dref[1] * drug_B / dref[2]) |
          group_id,
      data = codata_combo2,
      prior_EX_mu_comp = list(
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))
      ),
      prior_EX_tau_mean_comp = matrix(
        0,
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_tau_sd_comp = matrix(
        c(log(4) / 1.96, log(4) / 1.96, log(4) / 1.96, log(4) / 1.96),
        nrow = num_comp,
        ncol = 2,
        byrow = TRUE
      ),
      prior_EX_mu_mean_inter = 0,
      prior_EX_mu_sd_inter = 1.121,
      prior_EX_tau_mean_inter = matrix(0, nrow = num_strata, ncol = num_inter),
      prior_EX_tau_sd_inter = matrix(
        log(4) / 1.96,
        nrow = num_strata,
        ncol = num_inter
      ),
      prior_is_EXNEX_comp = rep(FALSE, num_comp),
      prior_is_EXNEX_inter = rep(FALSE, num_inter),
      prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = num_comp),
      prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = num_inter),
      ## setup sampling to be fast, but do not trigger any rstan warning
      ## as the rstan warning would otherwise
      prior_tau_dist = 0,
      prior_PD = TRUE,
      init = 0,
      chains = 1,
      iter = 1000,
      warmup = 500,
      control = list(adapt_delta = 0.85)
    )),
    regexp = "and not deprecated arguments"
  )

  expect_no_error(
    blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug_A / dref[1])) |
          1 + I(log(drug_B / dref[2])) |
          0 + I(drug_A / dref[1] * drug_B / dref[2]) |
          group_id,
      data = codata_combo2,
      prior_EX_mu_comp = list(
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))
      ),
      prior_EX_tau_comp = list(list(
        mixmvnorm(c(1, 0, 0, diag(c(log(4) / 1.96, log(4) / 1.96)^2))),
        mixmvnorm(c(1, 0, 0, diag(c(log(4) / 1.96, log(4) / 1.96)^2)))
      )),
      prior_EX_mu_inter = mixmvnorm(c(1, 0, 1.121^2)),
      prior_EX_tau_inter = replicate(
        num_strata,
        mixmvnorm(c(
          1,
          rep.int(0, num_inter),
          diag(rep.int((log(4) / 1.96)^2, num_inter), num_inter)
        )),
        FALSE
      ),
      prior_is_EXNEX_comp = rep(FALSE, num_comp),
      prior_is_EXNEX_inter = rep(FALSE, num_inter),
      prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = num_comp),
      prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = num_inter),
      ## setup sampling to be fast, but do not trigger any rstan warning
      ## as the rstan warning would otherwise
      prior_tau_dist = 0,
      prior_PD = TRUE,
      init = 0,
      chains = 1,
      iter = 1000,
      warmup = 500,
      control = list(adapt_delta = 0.85)
    )
  )
})

test_that("blrm_exnex exits gracefully with an error message when Stan does not sample with rstan", {
  skip_on_cran()
  single_agent <- load_single_agent()

  ## give wrong argument to Stan sampler...should give us a proper error message
  o <- capture.output(
    expect_error(
      suppressMessages(with(
        single_agent,
        backend = "rstan",
        update(blrmfit, control = list(aBapt_delta = 0.85))
      )),
      "Calling Stan failed."
    ),
    type = "message"
  )
})

test_that("blrm_exnex exits gracefully with an error message when Stan does not sample with cmdstanr", {
  skip_if_no_cmdstan()
  skip_on_cran()
  single_agent <- load_single_agent()

  ## give wrong argument to Stan sampler...should give us a proper error message
  expect_error(
    with(
      single_agent,
      update(blrmfit, backend = "cmdstanr", control = list(adapt_delta = -0.85))
    ),
    "Assertion on"
  )
})

test_that("blrm_exnex accepts mixture priors with differing number of components", {
  skip_on_cran()
  combo2 <- load_combo2()

  withr::local_seed(67894356)

  expect_no_error(
    mix_comp_fit <- with(
      combo2,
      update(
        blrmfit,
        prior_PD = TRUE,
        init = 0,
        iter = 2000,
        warmup = 1000,
        chains = 1,
        prior_EX_mu_comp = list(
          mixmvnorm(
            c(0.5, logit(0.2), 0, diag(c(2^2, 1))),
            c(0.5, logit(0.2), 0, diag(c(2^2, 1)))
          ),
          mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))
        )
      )
    )
  )
})

test_that("blrm_exnex does not accept hierarchical prior arguments whenever the hierarchical model is disabled.", {
  skip_on_cran()

  dref <- c(6, 960)

  num_comp <- 2 # two investigational drugs
  num_inter <- 1 # one drug-drug interaction needs to be modeled
  num_groups <- nlevels(codata_combo2$group_id) # no stratification needed
  num_strata <- 1 # no stratification needed

  expect_error(
    blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug_A / dref[1])) |
          1 + I(log(drug_B / dref[2])) |
          0 + I(drug_A / dref[1] * drug_B / dref[2]) |
          group_id,
      data = codata_combo2,
      prior_EX_mu_comp = list(
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))
      ),
      prior_EX_tau_comp = list(list(
        mixmvnorm(c(1, 0, 0, diag(c(log(4) / 1.96, log(4) / 1.96)^2))),
        mixmvnorm(c(1, 0, 0, diag(c(log(4) / 1.96, log(4) / 1.96)^2)))
      )),
      prior_EX_mu_inter = mixmvnorm(c(1, 0, 1.121^2)),
      prior_EX_tau_inter = replicate(
        num_strata,
        mixmvnorm(c(
          1,
          rep.int(0, num_inter),
          diag(rep.int((log(4) / 1.96)^2, num_inter), num_inter)
        )),
        FALSE
      ),
      prior_is_EXNEX_comp = rep(FALSE, num_comp),
      prior_is_EXNEX_inter = rep(FALSE, num_inter),
      prior_EX_prob_comp = matrix(1, nrow = num_groups, ncol = num_comp),
      prior_EX_prob_inter = matrix(1, nrow = num_groups, ncol = num_inter),
      ## setup sampling to be fast, but do not trigger any rstan warning
      ## as the rstan warning would otherwise
      prior_tau_dist = NULL,
      prior_PD = TRUE,
      init = 0,
      chains = 1,
      iter = 1000,
      warmup = 500,
      control = list(adapt_delta = 0.85)
    ),
    regex = "Hierarchical model structure disabled"
  )
})

test_that("blrm_exnex does require presence of a single group whenever the hierarchical model is disabled.", {
  skip_on_cran()

  dref <- c(6, 960)

  num_comp <- 2 # two investigational drugs
  num_inter <- 1 # one drug-drug interaction needs to be modeled
  num_groups <- nlevels(codata_combo2$group_id) # no stratification needed
  num_strata <- 1 # no stratification needed

  expect_message(
    blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug_A / dref[1])) |
          1 + I(log(drug_B / dref[2])) |
          0 + I(drug_A / dref[1] * drug_B / dref[2]) |
          group_id,
      data = codata_combo2,
      prior_EX_mu_comp = list(
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))
      ),
      prior_EX_mu_inter = mixmvnorm(c(1, 0, 1.121^2)),
      ## setup sampling to be fast, but do not trigger any rstan warning
      ## as the rstan warning would otherwise
      prior_tau_dist = NULL,
      prior_PD = TRUE,
      init = 0,
      chains = 1,
      iter = 1000,
      warmup = 500,
      control = list(adapt_delta = 0.85)
    ),
    regex = "found more than one group"
  )
})

test_that("blrm_exnex returns tau posteriors which are zero whenever the hierarchical model is disabled.", {
  skip_on_cran()

  dref <- c(6, 960)

  suppressWarnings(
    simple_fit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug_A / dref[1])) |
          1 + I(log(drug_B / dref[2])) |
          0 + I(drug_A / dref[1] * drug_B / dref[2]) |
          group_id,
      data = mutate(codata_combo2, group_id = "trial"),
      prior_EX_mu_comp = list(
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1)))),
        mixmvnorm(c(1, logit(0.2), 0, diag(c(2^2, 1))))
      ),
      prior_EX_mu_inter = mixmvnorm(c(1, 0, 1.121^2)),
      prior_tau_dist = NULL,
      prior_PD = FALSE,
      init = 0,
      chains = 1,
      iter = 1000,
      warmup = 500,
      control = list(adapt_delta = 0.85)
    )
  )

  expect_true(all(
    as_draws_matrix(simple_fit, variable = c("tau_log_beta", "tau_eta")) == 0
  ))
})

test_that("blrm_exnex returns tau posteriors which are zero whenever the hierarchical model is disabled for a case with strata.", {
  skip_on_cran()

  dref <- c(500, 500, 1000)
  num_comp <- 3
  num_inter <- choose(3, 2) + 1
  num_strata <- nlevels(hist_combo3$stratum_id)
  num_groups <- nlevels(hist_combo3$group_id)

  suppressWarnings(
    blrmfit <- blrm_exnex(
      cbind(num_toxicities, num_patients - num_toxicities) ~
        1 +
          I(log(drug_A / dref[1])) |
          1 + I(log(drug_B / dref[2])) |
          1 + I(log(drug_C / dref[3])) |
          0 +
            I(drug_A / dref[1] * drug_B / dref[2]) +
            I(drug_A / dref[1] * drug_C / dref[3]) +
            I(drug_B / dref[2] * drug_C / dref[3]) +
            I(drug_A / dref[1] * drug_B / dref[2] * drug_C / dref[3]) |
          stratum_id / group_id,
      data = hist_combo3,
      prior_EX_mu_comp = replicate(
        num_comp,
        mixmvnorm(c(1, logit(1 / 3), 0, diag(c(2^2, 1)))),
        FALSE
      ),
      prior_EX_mu_inter = mixmvnorm(c(
        1,
        rep.int(0, num_inter),
        diag((rep.int(sqrt(2) / 2, num_inter))^2)
      )),
      prior_tau_dist = NULL,
      prior_PD = TRUE,
      init = 0,
      chains = 1,
      iter = 1000,
      warmup = 500,
      control = list(adapt_delta = 0.85)
    )
  )

  expect_true(all(
    as_draws_matrix(blrmfit, variable = c("tau_log_beta", "tau_eta")) == 0
  ))
})


test_that("blrm_exnex outputs samples of MAP priors when requested", {
  skip_on_cran()
  combo2 <- load_combo2()

  suppressWarnings(
    mapfit <- with(
      combo2,
      update(
        blrmfit,
        data = mutate(
          codata_combo2,
          group_id = factor(
            group_id,
            levels = c(levels(codata_combo2$group_id), "new")
          )
        ),
        prior_EX_prob_comp = matrix(
          1,
          nrow = nlevels(codata_combo2$group_id) + 1,
          ncol = 2
        ),
        prior_EX_prob_inter = matrix(
          1,
          nrow = nlevels(codata_combo2$group_id) + 1,
          ncol = 1
        ),
        sample_map = TRUE
      )
    )
  )

  draws <- as_draws_rvars(mapfit)

  expect_choice("map_log_beta", variables(draws))
  expect_choice("map_eta", variables(draws))

  expect_equal(
    dimnames(draws$map_log_beta),
    list(
      NULL,
      c("I(log(drug_A/dref[1]))", "I(log(drug_B/dref[2]))"),
      c("intercept", "log_slope")
    )
  )
  expect_equal(
    dimnames(draws$map_eta),
    list(NULL, c("I(drug_A/dref[1] * drug_B/dref[2])"))
  )

  ## check that the sampled MAP corresponds about what we expect it to
  ## be
  alt_sample_map_comp <- draws$beta_group["new", , ]
  alt_sample_map_comp[,, "slope"] <- log(alt_sample_map_comp[,, "slope"])
  alt_sample_map_inter <- draws$eta_group["new", ]

  delta <- alt_sample_map_comp - draws$map_log_beta
  delta_sum <- summarise_draws(delta, "mean", "sd")
  expect_numeric(
    delta_sum$mean / delta_sum$sd,
    lower = qnorm(1E-2),
    upper = qnorm(1 - 1E-2),
    any.missing = FALSE,
    len = 4
  )
})
