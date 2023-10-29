set.seed(1)
n <- 1500
p <- 3
x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
delta <- rnorm(p)
theta_0 <- c(1, 1.2, -1)
y <- c(
  rpois(300, exp(x[1:300, ] %*% theta_0)),
  rpois(400, exp(x[301:700, ] %*% (theta_0 + delta))),
  rpois(300, exp(x[701:1000, ] %*% theta_0)),
  rpois(100, exp(x[1001:1100, ] %*% (theta_0 - delta))),
  rpois(200, exp(x[1101:1300, ] %*% theta_0)),
  rpois(200, exp(x[1301:n, ] %*% (theta_0 + delta)))
)
result <- fastcpd(
  formula = y ~ . - 1,
  data = data.frame(y = y, x = x),
  beta = (p + 1) * log(n) / 2,
  k = function(x) 0,
  family = "poisson",
  epsilon = 1e-5
)
summary(result)
