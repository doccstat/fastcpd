set.seed(1)
data <- c(
  rnorm(3000, 0, 1),
  rnorm(1000, 10, 1),
  rnorm(3000, 10, 20),
  rnorm(1000, 0, 1)
)
system.time(result <- fastcpd.mv(data))
summary(result)
plot(result)
