withr::local_seed(123105)
very_fast_sampling()

fixture <- withr::with_options(
  list(OncoBayes2.MC.save_warmup = TRUE),
  {
    example <- run_example("combo2")
    example$blrmfit
  }
)
