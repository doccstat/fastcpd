set.seed(1)
p <- 3
data <- if (requireNamespace("mvtnorm", quietly = TRUE)) {
  rbind(
    mvtnorm::rmvnorm(2e+5, mean = rep(0, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(1e+5, mean = rep(50, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(2e+5, mean = rep(0, p), sigma = diag(100, p)),
    mvtnorm::rmvnorm(2e+5, mean = rep(0, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(1e+5, mean = rep(50, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(2e+5, mean = rep(50, p), sigma = diag(100, p))
  )
} else {
  rbind(
    matrix(rnorm(p * 2e+5, mean = 0, sd = 1), ncol = p),
    matrix(rnorm(p * 1e+5, mean = 50, sd = 1), ncol = p),
    matrix(rnorm(p * 2e+5, mean = 0, sd = 10), ncol = p),
    matrix(rnorm(p * 2e+5, mean = 0, sd = 1), ncol = p),
    matrix(rnorm(p * 1e+5, mean = 50, sd = 1), ncol = p),
    matrix(rnorm(p * 2e+5, mean = 50, sd = 10), ncol = p)
  )
}
(result_time <- system.time(result <- fastcpd.mv(data)))
summary(result)
result@thetas[seq_len(p), ]
lapply(result@thetas[seq_len(p^2) + p, ], function(thetas) matrix(thetas, p))
