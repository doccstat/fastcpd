#' Solve poisson regression using Gradient Descent Extension to the
#' multivariate case
#'
#' @param data_new TODO
#' @param coef TODO
#' @param cum_coef TODO
#' @param cmatrix TODO
#' @param epsilon TODO
#' @param G TODO
#' @param L TODO
#' @param H TODO
#'
#' @return TODO
#' @export
cost_poisson_update <- function(data_new, coef, cum_coef, cmatrix, epsilon = 0.001, G = 10^10, L = -20, H = 20) {
  p <- length(data_new) - 1
  X_new <- data_new[1:p]
  Y_new <- data_new[p + 1]
  eta <- X_new %*% coef
  mu <- exp(eta)
  cmatrix <- cmatrix + (X_new %o% X_new) * min(as.numeric(mu), G)
  lik_dev <- as.numeric(-(Y_new - mu)) * X_new
  coef <- coef - solve(cmatrix + epsilon * diag(1, p), lik_dev)
  coef <- DescTools::Winsorize(coef, minval = L, maxval = H)
  cum_coef <- cum_coef + coef
  return(list(coef, cum_coef, cmatrix))
}
