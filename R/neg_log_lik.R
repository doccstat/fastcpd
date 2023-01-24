#' Solve logistic/poisson regression using Gradient Descent Extension to the
#' multivariate case
#'
#' @param data TODO
#' @param b TODO
#' @param family TODO
#' @param lambda TODO
#'
#' @return TODO
#' @export
neg_log_lik <- function(data, b, family, lambda = NULL) {
  p <- dim(data)[2] - 1
  x <- data[, 1:p, drop = FALSE]
  y <- data[, p + 1, drop = FALSE]
  xb <- x %*% b
  if (family == "binomial") {
    u <- as.numeric(xb)
    log_likelihood <- -y * u + log(1 + exp(u))
    return(sum(log_likelihood))
  } else if (family == "poisson") {
    u <- as.numeric(xb)
    log_likelihood <- -y * u + exp(u) + lfactorial(y)
    return(sum(log_likelihood))
  } else if (family == "gaussian") {
    stopifnot(!is.null(lambda))
    return(sum((y - xb)^2) / 2 + lambda * sum(abs(b)))
  }
}
