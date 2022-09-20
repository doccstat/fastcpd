#' Solve logistic regression using Gradient Descent Extension to the
#' multivariate case
#'
#' @param data_new TODO
#' @param coef TODO
#' @param cum_coef TODO
#' @param cmatrix TODO
#' @param epsilon TODO
#'
#' @return TODO
#' @export
cost_logistic_update <- function(data_new, coef, cum_coef,
                                 cmatrix, epsilon = 1e-10) {
  p <- length(data_new) - 1
  X_new <- data_new[1:p]
  Y_new <- data_new[p + 1]
  eta <- X_new %*% coef
  mu <- 1 / (1 + exp(-eta))
  cmatrix <- cmatrix + (X_new %o% X_new) * as.numeric((1 - mu) * mu)
  lik_dev <- as.numeric(-(Y_new - mu)) * X_new
  coef <- coef - solve(cmatrix + epsilon * diag(1, p), lik_dev)
  cum_coef <- cum_coef + coef
  return(list(coef, cum_coef, cmatrix))
}
