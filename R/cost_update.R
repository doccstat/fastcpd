#' Solve logistic regression using Gradient Descent Extension to the
#' multivariate case
#'
#' @param data TODO
#' @param theta_hat TODO
#' @param theta_sum TODO
#' @param hessian TODO
#' @param family TODO
#' @param momentum TODO
#' @param momentum_coef TODO
#' @param args_list TODO
#'
#' @return TODO
#' @export
cost_update <- function(
  data,
  theta_hat,
  theta_sum,
  hessian,
  family,
  momentum,
  momentum_coef,
  args_list
) {
  p <- length(data) - 1
  x <- data[-1]
  y <- data[1]
  if (family %in% c("binomial", "poisson")) {
    if (is.null(args_list$epsilon)) {
      epsilon <- 1e-10
    }
    if (family == "poisson") {
      if (is.null(args_list$G)) {
        G <- 10^10
      }
      if (is.null(args_list$L)) {
        L <- -20
      }
      if (is.null(args_list$H)) {
        H <- 20
      }
    }
  } else if (family == "gaussian") {
    stopifnot(!is.null(args_list$lambda))
  } else {
    stop("family must be one of binomial, poisson, gaussian")
  }

  if (family == "binomial") {
    prob <- 1 / (1 + exp(-x %*% theta_hat))
    hessian <- hessian + (x %o% x) * as.numeric((1 - prob) * prob)
    gradient <- as.numeric(-(y - prob)) * x
    momentum <- momentum_coef * momentum - solve(hessian + epsilon * diag(1, p), gradient)
    theta_hat <- theta_hat + momentum
  } else if (family == "poisson") {
    prob <- exp(x %*% theta_hat)
    hessian <- hessian + (x %o% x) * min(as.numeric(prob), G)
    gradient <- as.numeric(-(y - prob)) * x
    momentum <- momentum_coef * momentum - solve(hessian + epsilon * diag(1, p), gradient)
    theta_hat <- theta_hat + momentum
    theta_hat <- DescTools::Winsorize(theta_hat, minval = L, maxval = H)
  } else if (family == "gaussian") {
    prob <- x %*% theta_hat
    hessian <- hessian + x %o% x
    # B <- as.vector(cmatrix_inv%*%x)
    # cmatrix_inv <- cmatrix_inv - B%o%B/(1+sum(x*B))
    gradient <- as.numeric(-(y - prob)) * x
    momentum <- momentum_coef * momentum - solve(hessian, gradient)
    theta_hat <- theta_hat + momentum
    nc <- norm(hessian, type = "F") # the choice of norm affects the speed. Spectral norm is more accurate but slower than F norm.
    theta_hat <- sign(theta_hat) * pmax(abs(theta_hat) - args_list$lambda / nc, 0)
  } else {
    stop("family must be one of 'gaussian', 'binomial', or 'poisson'")
  }
  theta_sum <- theta_sum + theta_hat
  return(list(theta_hat, theta_sum, hessian, momentum))
}
