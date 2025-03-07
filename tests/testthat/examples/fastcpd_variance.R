if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(1)
  p <- 3
  data <- rbind(
    mvtnorm::rmvnorm(
      3e+5, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
    ),
    mvtnorm::rmvnorm(
      4e+5, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
    ),
    mvtnorm::rmvnorm(
      3e+5, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
    )
  )
  result_time <- system.time(
    result <- fastcpd.variance(data, r.progress = FALSE, cp_only = TRUE)
  )
  print(result_time)
  summary(result)
}
