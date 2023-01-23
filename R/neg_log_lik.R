#' Solve logistic/poisson regression using Gradient Descent Extension to the
#' multivariate case
#'
#' @param data TODO
#' @param b TODO
#'
#' @return TODO
#' @export
neg_log_lik <- function(data, b, family) {
  p <- dim(data)[2] - 1
  X <- data[, 1:p, drop = FALSE]
  Y <- data[, p + 1, drop = FALSE]
  u <- as.numeric(X %*% b)
  if (family == "binomial") {
    return(sum(-Y * u + log(1 + exp(u))))
  } else if (family == "poisson") {
    return(sum(-Y * u + exp(u) + lfactorial(Y)))
  }
}
