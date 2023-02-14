#' Generate data from logistic regression models with change-points
#'
#' @param n TODO
#' @param d TODO
#' @param tau TODO
#' @param Sigma TODO
#' @param family TODO
#' @param evar TODO
#'
#' @return TODO
#' @export
data_gen <- function(n, d, tau, Sigma, family, evar = NULL) {
  loc <- unique(c(0, tau, n))
  theta <- matrix(NA, nrow = d, ncol = length(tau) + 1)
  for (i in seq_len(length(tau) + 1)) {
    theta[, i] <- c((-1) ** (i + 1), rnorm(d - 1, 0.5, 0.1))
  }
  x <- mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = Sigma)
  y <- NULL
  for (i in 1:(length(tau) + 1)) {
    Xb <- x[(loc[i] + 1):loc[i + 1], , drop = FALSE] %*% theta[, i, drop = FALSE]
    y_new <- if (family == "binomial") {
      stats::rbinom(length(Xb), 1, 1 / (1 + exp(-Xb)))
    } else if (family == "poisson") {
      stats::rpois(length(Xb), exp(Xb))
    } else if (family == "gaussian") {
      stats::rnorm(length(Xb), sd = sqrt(evar)) + Xb
    } else {
      stop("family not supported")
    }
    y <- c(y, y_new)
  }
  data <- cbind(y, x)
  segment <- rep(seq_len(length(tau) + 1), diff(loc))
  return(list(data = data, theta = theta, segment = segment))
}
