set.seed(42)
n <- 400
x <- rep(1, n)
eps <- rnorm(n)
eps[c(50, 51, 52, 300, 301, 302)] <- -15
y <- c(0 + eps[1:200], 3 + eps[201:400])
dat <- cbind(y, x)
# lm is dominated by the six -15 outliers and reports spurious breaks
result_lm <- fastcpd_lm(dat, r.progress = FALSE)
# Quantile regression (tau=0.5) is robust to outliers below 50% and
# correctly finds only the true change at 200.
result_qr <- fastcpd_quantile(dat, order = 0.5, r.progress = FALSE)
