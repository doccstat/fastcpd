if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(1)
  p <- 3
  result <- fastcpd.mv(
    rbind(
      mvtnorm::rmvnorm(600, mean = rep(0, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(800, mean = rep(100, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(600, mean = rep(0, p), sigma = diag(900, p)),
      mvtnorm::rmvnorm(600, mean = rep(0, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(800, mean = rep(100, p), sigma = diag(1, p)),
      mvtnorm::rmvnorm(600, mean = rep(100, p), sigma = diag(900, p))
    )
  )
  summary(result)
}
