if (!requireNamespace("mvtnorm", quietly = TRUE)) utils::install.packages(
  "mvtnorm", repos = "https://cloud.r-project.org", quiet = TRUE
)

set.seed(1)
p <- 3
result <- fastcpd.variance(
  rbind(
    mvtnorm::rmvnorm(300, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))),
    mvtnorm::rmvnorm(400, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p))),
    mvtnorm::rmvnorm(300, rep(0, p), crossprod(matrix(runif(p^2) * 2 - 1, p)))
  )
)
summary(result)
