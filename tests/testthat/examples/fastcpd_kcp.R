set.seed(1)
p_values <- c(runif(200, 0, 0.05), runif(200, 0, 1))
result <- detect_kernel(p_values)
summary(result)

set.seed(7)
identical_result <- detect_kernel(
  rep(3, 20), order = c(8, 0), beta = 2, trim = 0, cp_only = TRUE
)
