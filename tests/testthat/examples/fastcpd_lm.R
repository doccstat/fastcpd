set.seed(1)
n <- 300
p <- 4
x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
theta_0 <- rbind(c(1, 3.2, -1, 0), c(-1, -0.5, 2.5, -2), c(0.8, 0, 1, 2))
y <- c(
  x[1:100, ] %*% theta_0[1, ] + rnorm(100, 0, 3),
  x[101:200, ] %*% theta_0[2, ] + rnorm(100, 0, 3),
  x[201:n, ] %*% theta_0[3, ] + rnorm(100, 0, 3)
)
result_lm <- fastcpd.lm(cbind(y, x))
summary(result_lm)
plot(result_lm)

set.seed(1)
n <- 600
p <- 4
d <- 2
x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
theta_1 <- matrix(runif(8, -3, -1), nrow = p)
theta_2 <- matrix(runif(8, -1, 3), nrow = p)
y <- rbind(
  x[1:350, ] %*% theta_1 + mvtnorm::rmvnorm(350, rep(0, d), 3 * diag(d)),
  x[351:n, ] %*% theta_2 + mvtnorm::rmvnorm(250, rep(0, d), 3 * diag(d))
)
result_mlm <- fastcpd.lm(cbind.data.frame(y = y, x = x), p.response = 2)
summary(result_mlm)
