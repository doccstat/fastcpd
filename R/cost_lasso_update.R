#' Solve logistic regression using Gradient Descent Extension to the
#' multivariate case
#'
#' @param data TODO
#' @param coef TODO
#' @param cum_coef TODO
#' @param cmatrix TODO
#' @param lambda TODO
#'
#' @return TODO
#' @export
cost_lasso_update <- function(data, coef, cum_coef, cmatrix, lambda) {
  p <- length(data) - 1
  x <- data[1:p]
  y <- data[p + 1]
  mu <- x %*% coef
  cmatrix <- cmatrix + x %o% x
  # B <- as.vector(cmatrix_inv%*%x)
  # cmatrix_inv <- cmatrix_inv - B%o%B/(1+sum(x*B))
  lik_dev <- as.numeric(-(y - mu)) * x
  coef <- coef - solve(cmatrix, lik_dev)
  nc <- norm(cmatrix, type = "F") # the choice of norm affects the speed. Spectral norm is more accurate but slower than F norm.
  coef <- sign(coef) * pmax(abs(coef) - lambda / nc, 0)
  cum_coef <- cum_coef + coef
  return(list(coef, cum_coef, cmatrix))
}
