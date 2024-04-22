if (requireNamespace("mvtnorm", quietly = TRUE)) {
  set.seed(1)
  p <- 3
  result <- fastcpd.variance(
    rbind(
      mvtnorm::rmvnorm(
        300, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
      ),
      mvtnorm::rmvnorm(
        400, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
      ),
      mvtnorm::rmvnorm(
        300, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))
      )
    )
  )
  summary(result)
}
