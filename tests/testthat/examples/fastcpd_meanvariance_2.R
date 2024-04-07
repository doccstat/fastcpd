set.seed(1)
p <- 4
result <- fastcpd.mv(
  rbind(
    mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(100, p)),
    mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(300, mean = rep(10, p), sigma = diag(100, p))
  )
)
summary(result)
