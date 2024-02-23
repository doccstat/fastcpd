for (package in c("ggplot2", "mvtnorm")) {
  if (!requireNamespace(package, quietly = TRUE)) utils::install.packages(
    package, repos = "https://cloud.r-project.org", quiet = TRUE
  )
}

set.seed(1)
p <- 4
x <- mvtnorm::rmvnorm(400, rep(0, p), diag(p))
theta <- rbind(rnorm(p, 0, 1), rnorm(p, 2, 1))
y <- c(
  rbinom(150, 1, 1 / (1 + exp(-x[1:150, ] %*% theta[1, ]))),
  rbinom(250, 1, 1 / (1 + exp(-x[151:400, ] %*% theta[2, ])))
)
result <- suppressWarnings(
  fastcpd.binomial(cbind(y, x), cost_adjustment = NULL)
)
summary(result)
plot(result)
