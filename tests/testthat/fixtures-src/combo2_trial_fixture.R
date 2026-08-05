withr::local_seed(123103)
very_fast_sampling()

fixture <- with(
  examples$combo2,
  suppressWarnings(suppressMessages(
    blrm_trial(
      histdata,
      dose_info,
      drug_info,
      simplified_prior = TRUE
    )
  ))
)
