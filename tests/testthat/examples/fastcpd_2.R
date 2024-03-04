if (!requireNamespace("mvtnorm", quietly = TRUE)) utils::install.packages(
  "mvtnorm", repos = "https://cloud.r-project.org", quiet = TRUE
)

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
result_mlm <- fastcpd(
  cbind(y.1, y.2) ~ . - 1, cbind.data.frame(y = y, x = x), family = "lm"
)
summary(result_mlm)
