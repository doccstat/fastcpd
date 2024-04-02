set.seed(1)
data <- c(rnorm(300, 0, 1), rnorm(400, 0, 100), rnorm(300, 0, 1))
result <- fastcpd.variance(data)
summary(result)
