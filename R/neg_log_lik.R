#' Solve logistic regression using Gradient Descent Extension to the
#' multivariate case
#'
#' @param data TODO
#' @param b TODO
#'
#' @return TODO
#' @export
neg_log_lik <- function(data, b) {
  p <- dim(data)[2] - 1
  X <- data[, 1:p, drop = FALSE]
  Y <- data[, p + 1, drop = FALSE]
  u <- as.numeric(X %*% b)
  L <- -Y * u + log(1 + exp(u))
  return(sum(L))
}
