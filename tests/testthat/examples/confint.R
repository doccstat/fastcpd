# The well-log data contain outliers, so a robust variance estimate keeps
# the detection and the confidence intervals stable.
sigma2 <- estimate_variance_median(well_log)
result <- detect_mean(well_log, trim = 0.001, variance_estimation = sigma2)

(cp_profile_interval <- confint(
  result,
  parm = "cp",
  method = "profile",
  level = 0.8,
  window = 8
))

(theta_wald_interval <- confint(
  result,
  parm = "theta",
  method = "wald"
))

(cp_bootstrap_interval <- confint(
  result,
  parm = "cp",
  method = "bootstrap",
  B = 5,
  seed = 10
))

variance_result <- detect_variance(well_log)
# Wald intervals are not available for variance-family fits.
