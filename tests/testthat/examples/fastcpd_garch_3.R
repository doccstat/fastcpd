set.seed(1)
n <- 300
sigma_2 <- rep(1, n + 1)
x <- rep(0, n + 1)
for (i in seq_len(100)) {
  sigma_2[i + 1] <- 10 + 0.5 * x[i]^2 + 0.3 * sigma_2[i]
  x[i + 1] <- rnorm(1, 0, sqrt(sigma_2[i + 1]))
}
for (i in 101:n) {
  sigma_2[i + 1] <- 0.2 + 0.05 * x[i]^2 + 0.1 * sigma_2[i]
  x[i + 1] <- rnorm(1, 0, sqrt(sigma_2[i + 1]))
}
result <- suppressWarnings(
  fastcpd.garch(x[-1], c(1, 1), include.mean = FALSE)
)
summary(result)
plot(result)
