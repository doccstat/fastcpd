set.seed(1)
p <- 1
result <- fastcpd.meanvariance(c(
  rnorm(300, 0, 1),
  rnorm(400, 10, 1),
  rnorm(300, 0, 10),
  rnorm(300, 0, 1),
  rnorm(400, 10, 1),
  rnorm(300, 10, 10)
))
summary(result)
plot(result)
