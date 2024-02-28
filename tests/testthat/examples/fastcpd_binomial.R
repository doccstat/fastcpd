for (package in c("ggplot2", "mvtnorm")) {
  if (!requireNamespace(package, quietly = TRUE)) utils::install.packages(
    package, repos = "https://cloud.r-project.org", quiet = TRUE
  )
}

set.seed(1)
n <- 500
p <- 4
x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
theta <- rbind(rnorm(p, 0, 1), rnorm(p, 2, 1))
y <- c(
  rbinom(300, 1, 1 / (1 + exp(-x[1:300, ] %*% theta[1, ]))),
  rbinom(200, 1, 1 / (1 + exp(-x[301:n, ] %*% theta[2, ])))
)
result <- suppressWarnings(fastcpd.binomial(cbind(y, x)))
summary(result)
plot(result)
