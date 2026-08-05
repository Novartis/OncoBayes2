withr::local_seed(123107)
very_fast_sampling()

example <- run_example("combo3")
suppressWarnings(
  fixture <- with(example, update(blrmfit, sample_map = TRUE))
)
