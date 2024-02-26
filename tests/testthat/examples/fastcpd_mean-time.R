# Test `fastcpd.mean()` running time. The running time is expected to be less
# than 1.5 seconds.
set.seed(1)
data <- c(rnorm(10000), rnorm(1000) + 1)
(result_time <- system.time(
  result <- fastcpd.mean(data, r.progress = FALSE, cp_only = TRUE)
))
result@cp_set
