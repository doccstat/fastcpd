cross_language_mean_data <- c(rep(0, 50), rep(5, 50))
cross_language_mean_result <- fastcpd.mean(
  cross_language_mean_data,
  beta = 5,
  cost_adjustment = "BIC",
  trim = 0,
  variance_estimation = 1
)
