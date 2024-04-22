if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(1)
  p <- 1
  x <- mvtnorm::rmvnorm(300, rep(0, p), diag(p))
  theta_0 <- matrix(c(1, -1, 0.5))
  y <- c(
    x[1:100, ] * theta_0[1, ] + rnorm(100, 0, 1),
    x[101:200, ] * theta_0[2, ] + rnorm(100, 0, 1),
    x[201:300, ] * theta_0[3, ] + rnorm(100, 0, 1)
  )
  result <- fastcpd.lm(cbind(y, x), r.clock = "fastcpd_profiler")
  summary(result)
  plot(result)

  if (requireNamespace("RcppClock", quietly = TRUE)) {
    library(RcppClock)
    plot(fastcpd_profiler)
  }
}
