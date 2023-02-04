```r
library(cpd)
set.seed(100)
n <- 1000
p <- 3
tau <- c(0.3, 0.85, 1)
mu <- matrix(c(-1, rnorm(p - 1, 0.5, 0.1)), p, n)
K <- length(tau)
for (i in 1:(K - 1))
{
  st <- floor(n * tau[i] + 1)
  en <- floor(n * tau[i + 1])
  mu[, st:en] <- c((-1)^(i + 1), rnorm(p - 1, 0.5, 0.1))
}
g_tr <- rep(1:K, diff(c(0, tau * n)))

# The choice of beta = (p+1)*log leads to the BIC

beta <- (p + 1) * log(n) / 2

reptim <- 10
rand_gd <- rand_va <- rep(NA, reptim)
for (i in 1:reptim) {
  X <- matrix(rnorm(n * (p - 1), 0), n, p - 1)
  X <- cbind(1, X)
  eta <- apply(t(mu) * X, 1, sum)
  prob <- 1 / (1 + exp(-eta))
  Y <- rbinom(n, size = 1, prob = prob)
  data <- cbind(Y, X)

  # seGD
  cp_set <- CP(data, beta, family = "binomial")$cp
  cp_gd <- cp_set[!(cp_set == 0)]
  K_est <- length(cp_gd) + 1
  cp_un <- unique(c(0, cp_gd, n))
  g_est <- rep(1:K_est, diff(cp_un))
  rand_gd[i] <- fossil::rand.index(g_tr, g_est)

  # vanilla
  cp_set <- CP_vanilla(data, beta, family = "binomial")$cp
  cp_va <- cp_set[!(cp_set == 0)]
  K_est <- length(cp_va) + 1
  cp_un <- unique(c(0, cp_va, n))
  g_est <- rep(1:K_est, diff(cp_un))
  rand_va[i] <- rand.index(g_tr, g_est)

  print(i)
}

c(mean(rand_gd), mean(rand_va))

n <- 1000
d <- 2
Sigma <- diag(1, d)
true.cp.loc <- c(200, 300, 400)
seg <- length(true.cp.loc) + 1
true.coef <- matrix(rnorm(seg * d), d, seg)
data_gen(n, d, true.coef, true.cp.loc, Sigma, "binomial")








n <- 1500
d <- 5
rho <- 0.9
Sigma <- array(0, c(d, d))
for (i in 1:d)
{
  Sigma[i, ] <- rho^(abs(i - (1:d)))
}

delta <- c(5, 7, 9, 11, 13)
a.sq <- 1
delta.new <- delta * sqrt(a.sq) / sqrt(as.numeric(t(delta) %*% Sigma %*% delta))
true.cp.loc <- c(375, 750, 1125)

true.coef <- matrix(0, nrow = d, ncol = length(true.cp.loc) + 1)
true.coef[, 1] <- c(1, 1.2, -1, 0.5, -2)
true.coef[, 2] <- true.coef[, 1] + delta.new
true.coef[, 3] <- true.coef[, 1]
true.coef[, 4] <- true.coef[, 3] - delta.new

out <- data_gen(n, d, true.coef, true.cp.loc, Sigma, "poisson")
data <- out[[1]]
g_tr <- out[[2]]
beta <- log(n) * (d + 1) / 2
CP(data, beta, trim = 0.03, B = 10, "poisson", epsilon = 0.001, G = 10^10, L = -20, H = 20)

CP_vanilla(data, beta, family = "poisson")









n <- 1000
s <- 3
d <- 50
evar <- 0.5
Sigma <- diag(1, d)
true.cp.loc <- c(100, 300, 500, 800, 900)
seg <- length(true.cp.loc) + 1
true.coef <- matrix(rnorm(seg * s), s, seg)
true.coef <- rbind(true.coef, matrix(0, d - s, seg))
out <- data_gen(n, d, true.coef, true.cp.loc, Sigma, family = "gaussian", evar)
data <- out[[1]]
beta <- log(n) / 2 # beta here has different meaning

start_time <- Sys.time()
CP_vanilla(data, beta, family = "gaussian", B = 10)
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
CP(data, beta, B = 10, trim = 0.025, family = "gaussian", epsilon = 1e-5)
end_time <- Sys.time()
end_time - start_time



####### BUGS #######
# fastglm::fastglm(matrix(0, ncol = 0, nrow = 1), 1, family = binomial())
```