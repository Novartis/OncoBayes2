withr::local_seed(123104)
very_fast_sampling()

example <- run_example("combo2_trial")
suppressWarnings(
  fixture <- with(example$combo2_trial, update(blrmfit, sample_map = TRUE))
)
