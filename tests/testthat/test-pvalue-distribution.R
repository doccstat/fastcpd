set.seed(1)
data <- c(
  runif(100, 0, 1),
  runif(100, 0, 0.025),
  runif(100, 0, 1),
  runif(100, 0, 0.025),
  runif(100, 0, 1)
)
plot(fastcpd::fastcpd.binomial(cbind(data <= 0.05, 1)))

rank_signal <- numeric(length(data))
for (t in seq_along(rank_signal)) {
  rank_signal[t] <- sum(data <= data[t]) - (length(data) + 1) / 2
}
rank_signal_empirical_covariance <- sum((rank_signal + 0.5)^2) / length(data)
rank_cost <- function(x) {
  -length(x) * mean(x)^2 / rank_signal_empirical_covariance
}
fastcpd::fastcpd(
  formula = ~ . - 1,
  data = data.frame(rank_signal),
  cost = rank_cost
)

gamma <- 0.1
cost_rbf <- function(x) {
  length(x) - sum(exp(- gamma * as.matrix(stats::dist(x))^2)) / length(x)
}
fastcpd::fastcpd(
  formula = ~ . - 1,
  data = data.frame(data),
  cost = cost_rbf,
  beta = 1
)
