#' Solve poisson regression using Gradient Descent Extension to the
#' multivariate case
#'
#' @param data TODO
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
cost_poisson_update <- function(data,
                                coef,
                                cum_coef,
                                cmatrix,
                                epsilon = 0.001,
                                G = 10^10,
                                L = -20,
                                H = 20) {
  p <- length(data) - 1
  x <- data[1:p]
  y <- data[p + 1]
  eta <- x %*% coef
  mu <- exp(eta)
  cmatrix <- cmatrix + (x %o% x) * min(as.numeric(mu), G)
  lik_dev <- as.numeric(-(y - mu)) * x
  coef <- coef - solve(cmatrix + epsilon * diag(1, p), lik_dev)
  coef <- DescTools::Winsorize(coef, minval = L, maxval = H)
  cum_coef <- cum_coef + coef
  return(list(coef, cum_coef, cmatrix))
}
