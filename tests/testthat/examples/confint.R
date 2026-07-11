set.seed(1)
data <- c(stats::rnorm(40, 0, 0.2), stats::rnorm(40, 3, 0.2))
result <- fastcpd.mean(data)

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

set.seed(3)
variance_data <- c(stats::rnorm(30, 0, 1), stats::rnorm(30, 0, 3))
variance_result <- fastcpd.variance(variance_data)
# Wald intervals are not available for variance-family fits.
