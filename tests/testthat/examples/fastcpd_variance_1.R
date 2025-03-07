set.seed(1)
data <- c(rnorm(3e+5, 0, 1), rnorm(4e+5, 0, 5), rnorm(3e+5, 0, 1))
(result_time <- system.time(
  result <- fastcpd.variance(data, r.progress = FALSE, cp_only = TRUE)
))
result@cp_set
