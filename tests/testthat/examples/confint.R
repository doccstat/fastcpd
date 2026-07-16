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

# ARIMA profile intervals use the same segment-local native likelihood as the
# detector rather than a separate stats::arima refit.
arima_small <- rep(c(0.1, -0.1), 20)
arima_large <- rep(c(2, -2), length.out = 41)
arima_data <- c(0, cumsum(c(arima_small, arima_large)))
arima_result <- detect_arima(
  arima_data,
  order = c(0, 1, 0)
)
arima_profile_interval <- confint(
  arima_result,
  parm = "cp",
  method = "profile",
  level = 0.8,
  window = 1
)
