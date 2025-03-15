set.seed(1)
p <- 3
data <- rbind(
  matrix(rnorm(p * 3e+5, mean = 0, sd = 10), ncol = p),
  matrix(rnorm(p * 4e+5, mean = 50, sd = 10), ncol = p),
  matrix(rnorm(p * 3e+5, mean = 2, sd = 10), ncol = p)
)
system.time(result <- fastcpd.mean(data, r.progress = FALSE, cp_only = TRUE))
summary(result)
