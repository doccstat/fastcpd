if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(1)
  p <- 3
  data <- rbind(
    mvtnorm::rmvnorm(3e+5, mean = rep(0, p), sigma = diag(100, p)),
    mvtnorm::rmvnorm(4e+5, mean = rep(50, p), sigma = diag(100, p)),
    mvtnorm::rmvnorm(3e+5, mean = rep(2, p), sigma = diag(100, p))
  )
  result_time <- system.time(
    result <- fastcpd.mean(data, r.progress = FALSE, cp_only = TRUE)
  )
  print(result_time)
  summary(result)
}
