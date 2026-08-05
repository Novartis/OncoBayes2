withr::local_seed(123110)
very_fast_sampling()

fixture <- with(
  examples$combo2,
  {
    dose_info <- dplyr::mutate(dose_info, drug1 = 1.0 * drug1)
    suppressWarnings(
      blrm_trial(
        histdata,
        dose_info,
        drug_info,
        simplified_prior = TRUE,
        interval_prob = c(0, 0.16, 0.33, 1),
        interval_max_mass = c(under = 0.2, target = 1, over = 0.25)
      )
    )
  }
)
