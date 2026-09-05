numeric_beta_lm_data <- cbind(
  y = rep(c(0, 1), 5),
  x = rep(1, 10)
)

numeric_beta_lm_result <- detect_lm(
  numeric_beta_lm_data,
  beta = 1e6,
  cost_adjustment = "BIC",
  segment_count = 1,
  cp_only = TRUE
)

numeric_beta_lm_explicit_variance <- detect_lm(
  numeric_beta_lm_data,
  beta = 1e6,
  cost_adjustment = "BIC",
  segment_count = 1,
  variance_estimation = 1,
  cp_only = TRUE
)
