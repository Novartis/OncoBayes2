context("plot curve tests")

plot_trial_fixture_names <- c(
  single_agent = "single_agent_trial",
  combo2 = "combo2_trial",
  combo3 = "combo3_trial"
)

# plot_toxicity_curve.blrmfit ----------------------------------------------
test_that("plot_toxicity_curve.blrmfit works for single-agent example", {
  skip_on_cran()

  single_agent <- load_test_fixture("single_agent_example")

  expect_gg(plot_toxicity_curve(single_agent$blrmfit, x = vars(drug_A)))

  expect_gg(plot_toxicity_curve(
    single_agent$blrmfit,
    x = vars(drug_A),
    transform = FALSE
  ))

  expect_gg(plot_toxicity_curve(
    single_agent$blrmfit,
    x = vars(drug_A),
    transform = FALSE,
    ylim = c(0, 0.5)
  ))

  expect_gg(
    plot_toxicity_curve(
      single_agent$blrmfit,
      x = vars(drug_A),
      transform = FALSE,
      ylim = c(0, 0.5),
      xlim = c(0.1, 25)
    ) +
      scale_x_log10()
  )

  expect_gg(plot_toxicity_curve(
    single_agent$blrmfit,
    x = vars(drug_A),
    transform = FALSE,
    ylim = c(0, 0.5),
    xlim = c(0.1, 25)
  ))

  # check nonsense variable
  expect_error(plot_toxicity_curve(single_agent$blrmfit, x = vars(drug_B)))

  # check character
  expect_gg(plot_toxicity_curve(single_agent$blrmfit, x = "drug_A"))

  # check nonsense character
  expect_error(plot_toxicity_curve(single_agent$blrmfit, x = "drug_B"))

  # check passing two vars()
  expect_gg(plot_toxicity_curve(
    single_agent$blrmfit,
    x = vars(drug_A),
    group = vars(group_id)
  ))

  # one var(), one formula
  expect_gg(plot_toxicity_curve(
    single_agent$blrmfit,
    x = vars(drug_A),
    group = ~group_id
  ))

  # one var(), one nonsense formula
  expect_error(plot_toxicity_curve(
    single_agent$blrmfit,
    x = vars(drug_A),
    group = ~drug_B
  ))

  # duplicated arguments
  expect_error(plot_toxicity_curve(
    single_agent$blrmfit,
    x = vars(drug_A),
    group = vars(drug_A)
  ))

  # duplicated arguments
  expect_error(plot_toxicity_curve(
    single_agent$blrmfit,
    x = vars(drug_A),
    group = ~drug_A
  ))

  # newdata
  expect_gg(plot_toxicity_curve(
    single_agent$blrmfit,
    x = vars(drug_A),
    group = ~group_id,
    newdata = tibble(
      group_id = unique(single_agent$blrmfit$data$group_id),
      drug_A = c(0.1, 0.3)
    )
  ))
})


test_that("plot_toxicity_curve.blrmfit works for combo2 example", {
  skip_on_cran()

  combo2 <- load_test_fixture("combo2_example")

  expect_gg(plot_toxicity_curve(
    combo2$blrmfit,
    x = vars(drug_A),
    group = vars(group_id, drug_B)
  ))

  expect_gg(plot_toxicity_curve(
    combo2$blrmfit,
    x = "drug_A",
    group = c("group_id", "drug_B")
  ))

  expect_gg(plot_toxicity_curve(
    combo2$blrmfit,
    x = "drug_A",
    group = ~ group_id + drug_B
  ))

  expect_gg(plot_toxicity_curve(
    combo2$blrmfit,
    x = "drug_A",
    group = ~ group_id + drug_B,
    transform = FALSE
  ))

  expect_gg(plot_toxicity_curve(
    combo2$blrmfit,
    x = "drug_A",
    group = ~ group_id + drug_B,
    transform = FALSE,
    ylim = c(0, 0.5),
    xlim = c(0.1, 25)
  ))

  expect_gg(plot_toxicity_curve(
    combo2$blrmfit,
    x = "drug_A",
    group = ~ group_id + drug_B,
    transform = FALSE,
    newdata = expand_grid(
      group_id = factor("trial_AB", levels(combo2$blrmfit$data$group_id)),
      drug_A = c(1, 5),
      drug_B = c(200, 300, 400)
    )
  ))
})


test_that("plot_toxicity_curve.blrm_trial works for blrm_trial examples", {
  skip_on_cran()

  for (fixture_name in plot_trial_fixture_names) {
    trial <- load_test_fixture(fixture_name)

    expect_gg(plot_toxicity_curve(trial))

    expect_gg(plot_toxicity_curve(trial, ewoc_shading = FALSE))
  }
})


# Visual regression testing ----------------------------------------------------

test_that("plot_toxicity_curve.blrmfit renders single-agent plots correctly", {
  testthat::skip_on_cran()

  testthat::skip_if_not_installed("vdiffr", minimum_version = min_vdiffr)
  testthat::skip_if_not(identical(Sys.getenv("TEST_VDIFFR"), "true"))

  single_agent <- load_test_fixture("single_agent_example")

  # basic
  p <- plot_toxicity_curve(single_agent$blrmfit, x = vars(drug_A))
  vdiffr::expect_doppelganger("single-agent blrmfit", p)

  # change limits
  p <- plot_toxicity_curve(
    single_agent$blrmfit,
    x = vars(drug_A),
    xlim = c(10, 25),
    ylim = c(0.1, 0.5)
  )
  vdiffr::expect_doppelganger("single-agent blrmfit limits", p)

  # change width
  p <- plot_toxicity_curve(
    single_agent$blrmfit,
    x = vars(drug_A),
    prob = 0.25,
    prob_outer = 0.5
  )
  vdiffr::expect_doppelganger("single-agent blrmfit width", p)
})


test_that("plot_toxicity_curve.blrmfit renders combo2 plots correctly", {
  testthat::skip_on_cran()

  testthat::skip_if_not_installed("vdiffr", minimum_version = min_vdiffr)
  testthat::skip_if_not(identical(Sys.getenv("TEST_VDIFFR"), "true"))

  combo2 <- load_test_fixture("combo2_example")

  nd <- filter(
    select(combo2$blrmfit$data, -starts_with("num")),
    group_id == "trial_AB"
  )

  # basic
  p <- plot_toxicity_curve(
    combo2$blrmfit,
    x = vars(drug_A),
    group = ~ drug_B + group_id,
    newdata = nd
  )
  vdiffr::expect_doppelganger("combo2 blrmfit", p)

  # change variable mapping
  p <- plot_toxicity_curve(
    combo2$blrmfit,
    x = vars(drug_B),
    group = ~ drug_A + group_id,
    newdata = nd
  )
  vdiffr::expect_doppelganger("combo2 blrmfit variables", p)

  # add another group_id
  p <- plot_toxicity_curve(
    combo2$blrmfit,
    x = vars(drug_B),
    group = ~ drug_A + group_id,
    newdata = expand_grid(
      select(nd, -group_id),
      group_id = factor(
        c("trial_AB", "IIT"),
        levels(combo2$blrmfit$group_fct)
      )
    )
  )
  vdiffr::expect_doppelganger("combo2 blrmfit newdata", p)

  # change limits
  p <- plot_toxicity_curve(
    combo2$blrmfit,
    x = vars(drug_A),
    group = ~ drug_B + group_id,
    newdata = nd,
    xlim = c(5, 10),
    ylim = c(0.1, 0.5)
  )
  vdiffr::expect_doppelganger("combo2 blrmfit limits", p)

  # change width
  p <- plot_toxicity_curve(
    combo2$blrmfit,
    x = vars(drug_A),
    group = ~ drug_B + group_id,
    newdata = nd,
    prob = 0.25,
    prob_outer = 0.5
  )
  vdiffr::expect_doppelganger("combo2 blrmfit width", p)
})


for (ex in names(plot_trial_fixture_names)) {
  test_that(
    paste("plot_toxicity_curve.blrm_trial renders", ex, "plots correctly"),
    {
      testthat::skip_on_cran()

      testthat::skip_if_not_installed("vdiffr", minimum_version = min_vdiffr)
      testthat::skip_if_not(identical(Sys.getenv("TEST_VDIFFR"), "true"))

      example <- load_test_fixture(plot_trial_fixture_names[[ex]])

      # blrm_trial
      p <- plot_toxicity_curve(example)
      vdiffr::expect_doppelganger(paste(ex, "blrm_trial"), p)

      p <- plot_toxicity_curve(example, ewoc_shading = FALSE)
      vdiffr::expect_doppelganger(paste(ex, "blrm_trial no ewoc"), p)
    }
  )
}
