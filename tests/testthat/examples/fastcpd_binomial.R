if (!requireNamespace("ggplot2", quietly = TRUE)) utils::install.packages(
  "ggplot2", repos = "https://cloud.r-project.org", quiet = TRUE
)

set.seed(1)
x <- matrix(rnorm(1500, 0, 1), ncol = 5)
theta <- rbind(rnorm(5, 0, 1), rnorm(5, 2, 1))
y <- c(
  rbinom(125, 1, 1 / (1 + exp(-x[1:125, ] %*% theta[1, ]))),
  rbinom(175, 1, 1 / (1 + exp(-x[126:300, ] %*% theta[2, ])))
)
result <- suppressWarnings(fastcpd.binomial(cbind(y, x)))
summary(result)
plot(result)
