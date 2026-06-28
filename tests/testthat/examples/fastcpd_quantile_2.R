set.seed(42)
n <- 400
x <- rep(1, n)
eps <- rnorm(n)
eps[c(50, 51, 52, 300, 301, 302)] <- -15
y <- c(0 + eps[1:200], 3 + eps[201:400])
dat <- cbind(y, x)
result_lm <- fastcpd_lm(dat)
summary(result_lm)
result_qr <- fastcpd_quantile(dat, order = 0.5)
summary(result_qr)
