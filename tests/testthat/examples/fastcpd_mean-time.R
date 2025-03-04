set.seed(1)
n <- 10^7
data <- c(rnorm(n / 2), rnorm(n / 2, 50))
(result_time <- system.time(
  result <- fastcpd.mean(data, r.progress = FALSE, cp_only = TRUE, variance_estimation = 1)
))
result@cp_set
