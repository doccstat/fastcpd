#' Solve logistic regression using Gradient Descent Extension to the
#' multivariate case
#'
#' @param data TODO
#' @param theta_hat TODO
#' @param theta_sum TODO
#' @param hessian TODO
#' @param tau TODO
#' @param i TODO
#' @param k TODO
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
  tau,
  i,
  k,
  family,
  momentum,
  momentum_coef,
  args_list
) {
  p <- ncol(data) - 1
  new_data_x <- data[nrow(data), -1]
  new_data_y <- data[nrow(data), 1]

  if (family %in% c("binomial", "poisson")) {
    epsilon <- if (is.null(args_list$epsilon)) 1e-10 else args_list$epsilon
    if (family == "poisson") {
      G <- if (is.null(args_list$G)) 10^10 else args_list$G
      L <- if (is.null(args_list$L)) -20 else args_list$L
      H <- if (is.null(args_list$H)) 20 else args_list$H
    }
  } else if (family == "gaussian") {
    stopifnot(!is.null(args_list$lambda))
  } else {
    stop("family must be one of binomial, poisson, gaussian")
  }

  if (family == "binomial") {
    prob <- 1 / (1 + exp(-new_data_x %*% theta_hat[, i]))
    hessian[, , i] <- hessian[, , i] + (new_data_x %o% new_data_x) * as.numeric((1 - prob) * prob)
    gradient <- as.numeric(-(new_data_y - prob)) * new_data_x
    momentum <- momentum_coef * momentum - solve(hessian[, , i] + epsilon * diag(1, p), gradient)
    theta_hat[, i] <- theta_hat[, i] + momentum
  } else if (family == "poisson") {
    prob <- exp(new_data_x %*% theta_hat[, i])
    hessian[, , i] <- hessian[, , i] + (new_data_x %o% new_data_x) * min(as.numeric(prob), G)
    gradient <- as.numeric(-(new_data_y - prob)) * new_data_x
    momentum <- momentum_coef * momentum - solve(hessian[, , i] + epsilon * diag(1, p), gradient)
    theta_hat[, i] <- theta_hat[, i] + momentum
    theta_hat[, i] <- DescTools::Winsorize(theta_hat[, i], minval = L, maxval = H)
  } else if (family == "gaussian") {
    hessian[, , i] <- hessian[, , i] + new_data_x %o% new_data_x
    # B <- as.vector(cmatrix_inv%*%new_data_x)
    # cmatrix_inv <- cmatrix_inv - B%o%B/(1+sum(new_data_x*B))
    gradient <- as.numeric(-(new_data_y - new_data_x %*% theta_hat[, i])) * new_data_x
    momentum <- momentum_coef * momentum - solve(hessian[, , i], gradient)
    theta_hat[, i] <- theta_hat[, i] + momentum
    nc <- norm(hessian[, , i], type = "F") # the choice of norm affects the speed. Spectral norm is more accurate but slower than F norm.
    theta_hat[, i] <- sign(theta_hat[, i]) * pmax(abs(theta_hat[, i]) - args_list$lambda / nc, 0)
  } else {
    stop("family must be one of 'gaussian', 'binomial', or 'poisson'")
  }

  for (kk in 1 + seq_len(k - 1)) {
    for (j in (tau + 1):nrow(data)) {

      if (family == "binomial") {
        prob <- 1 / (1 + exp(-data[j, -1] %*% theta_hat[, i]))
        hessian[, , i] <- hessian[, , i] + (data[j, -1] %o% data[j, -1]) * as.numeric((1 - prob) * prob)
        gradient <- as.numeric(-(data[j, 1] - prob)) * data[j, -1]
        momentum <- momentum_coef * momentum - solve(hessian[, , i] + epsilon * diag(1, p), gradient)
        theta_hat[, i] <- theta_hat[, i] + momentum
      } else if (family == "poisson") {
        prob <- exp(data[j, -1] %*% theta_hat[, i])
        hessian[, , i] <- hessian[, , i] + (data[j, -1] %o% data[j, -1]) * min(as.numeric(prob), G)
        gradient <- as.numeric(-(data[j, 1] - prob)) * data[j, -1]
        momentum <- momentum_coef * momentum - solve(hessian[, , i] + epsilon * diag(1, p), gradient)
        theta_hat[, i] <- theta_hat[, i] + momentum
        theta_hat[, i] <- DescTools::Winsorize(theta_hat[, i], minval = L, maxval = H)
      } else if (family == "gaussian") {
        hessian[, , i] <- hessian[, , i] + data[j, -1] %o% data[j, -1]
        gradient <- as.numeric(-(data[j, 1] - data[j, -1] %*% theta_hat[, i])) * data[j, -1]
        momentum <- momentum_coef * momentum - solve(hessian[, , i], gradient)
        theta_hat[, i] <- theta_hat[, i] + momentum
        nc <- norm(hessian[, , i], type = "F") # the choice of norm affects the speed. Spectral norm is more accurate but slower than F norm.
        theta_hat[, i] <- sign(theta_hat[, i]) * pmax(abs(theta_hat[, i]) - args_list$lambda / nc, 0)
      } else {
        stop("family must be one of 'gaussian', 'binomial', or 'poisson'")
      }
    }
  }

  theta_sum[, i] <- theta_sum[, i] + theta_hat[, i]
  return(list(theta_hat[, i], theta_sum[, i], hessian[, , i], momentum))
}
