if (!requireNamespace("mvtnorm", quietly = TRUE)) utils::install.packages(
  "mvtnorm", repos = "https://cloud.r-project.org", quiet = TRUE
)

set.seed(1)
p <- 3
result <- fastcpd.mv(
  rbind(
    mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(50, p)),
    mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
    mvtnorm::rmvnorm(300, mean = rep(10, p), sigma = diag(50, p))
  )
)
summary(result)
