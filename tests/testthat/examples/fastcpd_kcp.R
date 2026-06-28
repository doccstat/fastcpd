set.seed(1)
p_values <- c(runif(200, 0, 0.05), runif(200, 0, 1))
result <- fastcpd_kcp(p_values)
summary(result)
