set.seed(1)
data <- matrix(c(
  rexp(300, rate = 1),
  rexp(400, rate = 0.1),
  rexp(300, rate = 1)
))
system.time(result <- detect_exponential(data))
summary(result)
plot(result)
