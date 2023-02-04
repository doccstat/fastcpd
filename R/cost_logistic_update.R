#' Solve logistic regression using Gradient Descent Extension to the
#' multivariate case
#'
#' @param data TODO
#' @param coef TODO
#' @param cum_coef TODO
#' @param cmatrix TODO
#' @param epsilon TODO
#'
#' @return TODO
#' @export
cost_logistic_update <- function(data, coef, cum_coef,
                                 cmatrix, epsilon = 1e-10) {
  p <- length(data) - 1
  x <- data[1:p]
  y <- data[p + 1]
  eta <- x %*% coef
  mu <- 1 / (1 + exp(-eta))
  cmatrix <- cmatrix + (x %o% x) * as.numeric((1 - mu) * mu)
  lik_dev <- as.numeric(-(y - mu)) * x
  coef <- coef - solve(cmatrix + epsilon * diag(1, p), lik_dev)
  cum_coef <- cum_coef + coef
  return(list(coef, cum_coef, cmatrix))
}
