if (!requireNamespace("mvtnorm", quietly = TRUE)) utils::install.packages(
  "mvtnorm", repos = "https://cloud.r-project.org", quiet = TRUE
)

set.seed(1)
p <- 3
data <- rbind(
  mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(100, p)),
  mvtnorm::rmvnorm(400, mean = rep(50, p), sigma = diag(100, p)),
  mvtnorm::rmvnorm(300, mean = rep(2, p), sigma = diag(100, p))
)
result <- fastcpd.mean(data)
summary(result)
