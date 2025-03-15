set.seed(1)
data <- matrix(c(
  rnorm(300, mean = 0, sd = 10),
  rnorm(400, mean = 50, sd = 10),
  rnorm(300, mean = 2, sd = 10)
))
system.time(result <- fastcpd.mean(data))
summary(result)
plot(result)
