withr::local_seed(56363)

dose <- c(5, 10, 20, 40, 60, 80)

drug_info <- tibble(
  drug_name = "Dose",
  dose_ref = 40,
  dose_unit = "unitless"
)

dose_info <- tibble(group_id = "trial", Dose = dose)

fixture <- suppressWarnings(suppressMessages(blrm_trial(
  data = NULL,
  drug_info = drug_info,
  dose_info = dose_info,
  prior_EX_mu_comp = mixmvnorm(c(1, logit(0.2), 0, diag(c(1, 0.7)^2))),
  prior_EX_tau_comp = mixmvnorm(c(1, 0, 0, diag(c(1, 1)^2))),
  prior_EX_prob_comp = matrix(1, nrow = 1, ncol = 1),
  prior_tau_dist = 0,
  prior_PD = FALSE,
  chains = 4,
  cores = 1,
  warmup = 1000,
  iter = 2000
)))
