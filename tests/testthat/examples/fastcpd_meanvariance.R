if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(1)
  p <- 3
  result_time <- system.time(result <- fastcpd.mv(
    rbind(
      mvtnorm::rmvnorm(2e+5, mean = rep(0, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(1e+5, mean = rep(50, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(2e+5, mean = rep(0, p), sigma = diag(100, p)),
      mvtnorm::rmvnorm(2e+5, mean = rep(0, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(1e+5, mean = rep(50, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(2e+5, mean = rep(50, p), sigma = diag(100, p))
    )
  ))
  print(result_time)
  summary(result)
}
