if (!requireNamespace("mvtnorm", quietly = TRUE)) utils::install.packages(
  "mvtnorm", repos = "https://cloud.r-project.org", quiet = TRUE
)

set.seed(1)
n <- 300
p <- 4
x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
theta <- rbind(c(1, 3.2, -1, 0), c(-1, -0.5, 2.5, -2), c(0.8, 0, 1, 2))
y <- c(
  x[1:100, ] %*% theta[1, ] + rnorm(100, 0, 3),
  x[101:200, ] %*% theta[2, ] + rnorm(100, 0, 3),
  x[201:300, ] %*% theta[3, ] + rnorm(100, 0, 3)
)
(sigma2 <- variance.lm(cbind(y, x)))

set.seed(1)
n <- 300
p <- 4
d <- 2
x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
theta <- cbind(c(1, 3.2, -1, 0), c(-1, -0.5, 2.5, -2), c(0.8, 0, 1, 2))
theta <- cbind(theta, theta)
y <- rbind(
  x[1:100, ] %*% theta[, 1:2] + mvtnorm::rmvnorm(100, rep(0, d), 3 * diag(d)),
  x[101:200, ] %*% theta[, 3:4] + mvtnorm::rmvnorm(100, rep(0, d), 3 * diag(d)),
  x[201:300, ] %*% theta[, 5:6] + mvtnorm::rmvnorm(100, rep(0, d), 3 * diag(d))
)
(sigma <- variance.lm(cbind(y, x), d = d))
