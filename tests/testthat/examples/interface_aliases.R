set.seed(20260712)
interface_data <- c(
  stats::rnorm(60, mean = 0, sd = 0.2),
  stats::rnorm(60, mean = 4, sd = 0.2)
)

detect_interface_result <- detect_mean(interface_data, beta = 5)
fastcpd_alias_result <- fastcpd_mean(interface_data, beta = 5)
detect_generic_result <- detect(
  formula = ~ x - 1,
  data = data.frame(x = interface_data),
  family = "mean",
  beta = 5
)

estimate_variance_interface_result <- estimate_variance_mean(interface_data)
variance_alias_result <- variance_mean(interface_data)
estimate_variance_generic_result <- estimate_variance(
  interface_data,
  family = "mean"
)
