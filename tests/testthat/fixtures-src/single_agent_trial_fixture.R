withr::local_seed(123112)
very_fast_sampling()

fixture <- with(
  examples$single_agent,
  suppressWarnings(suppressMessages(
    blrm_trial(
      histdata,
      dose_info,
      drug_info,
      simplified_prior = TRUE
    )
  ))
)
