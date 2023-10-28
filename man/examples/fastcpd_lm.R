set.seed(1)
p <- 4
n <- 300
cp <- c(100, 200)
x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
theta_0 <- rbind(c(1, 3.2, -1, 0), c(-1, -0.5, 2.5, -2), c(0.8, -0.3, 1, 1))
y <- c(
  x[1:cp[1], ] %*% theta_0[1, ] + rnorm(cp[1], 0, sd = 3),
  x[(cp[1] + 1):cp[2], ] %*% theta_0[2, ] + rnorm(cp[2] - cp[1], 0, sd = 3),
  x[(cp[2] + 1):n, ] %*% theta_0[3, ] + rnorm(n - cp[2], 0, sd = 3)
)
result_internal_variance <- fastcpd.lm(cbind(y, x))
