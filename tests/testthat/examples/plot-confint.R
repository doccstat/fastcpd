sigma2 <- estimate_variance_median(well_log)
result <- detect_mean(well_log, trim = 0.001, variance_estimation = sigma2)

cp_profile_interval <- confint(
  result,
  parm = "cp",
  method = "profile",
  level = 0.8,
  window = 8
)
plot(cp_profile_interval)

theta_wald_interval <- confint(result, parm = "theta", method = "wald")
plot(theta_wald_interval)
