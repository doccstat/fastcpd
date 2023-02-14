#' Update the cost values for the segmentation.
#'
#' @param data A data frame containing the data to be segmented.
#' @param theta_hat Estimated theta from the previous iteration.
#' @param theta_sum Sum of estimated theta from the previous iteration.
#' @param hessian Hessian matrix from the previous iteration.
#' @param tau Start of the current segment.
#' @param i Index of the current data in the whole data set.
#' @param k Number of epochs in SGD.
#' @param family Family of the model.
#' @param momentum Momentum from the previous iteration.
#' @param momentum_coef Momentum coefficient to be applied to the current
#'   momentum.
#' @param args_list Arguments to be passed to the model.
#' @param cost_gradient Gradient for custom cost function.
#' @param cost_hessian Hessian for custom cost function.
#'
#' @return A list containing new values of \code{theta_hat}, \code{theta_sum},
#'   \code{hessian}, and \code{momentum}.
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
    args_list,
    cost_gradient,
    cost_hessian) {
  p <- ncol(data) - 1
  new_data_x <- data[nrow(data), -1]
  new_data_y <- data[nrow(data), 1]

  if (family %in% c("binomial", "poisson")) {
    epsilon <- if (is.null(args_list$epsilon)) 1e-10 else args_list$epsilon
    if (family == "poisson") {
      G <- if (is.null(args_list$G)) 10^10 else args_list$G
      winsorise_minval <- if (is.null(args_list$L)) -20 else args_list$L
      winsorise_maxval <- if (is.null(args_list$H)) 20 else args_list$H
    }
  } else if (family == "gaussian") {
    stopifnot(!is.null(args_list$lambda))
    lambda <- args_list$lambda
  } else {
    stop("family must be one of binomial, poisson, gaussian")
  }

  if (family == "binomial") {
    hessian[, , i] <- cost_hessian(
      data[nrow(data), ], theta_hat[, i], hessian[, , i], family, NULL
    )
    gradient <- cost_gradient(data[nrow(data), ], theta_hat[, i], family)
    momentum_step <- solve(hessian[, , i] + epsilon * diag(1, p), gradient)
    momentum <- momentum_coef * momentum - momentum_step
    theta_hat[, i] <- theta_hat[, i] + momentum
  } else if (family == "poisson") {
    hessian[, , i] <- cost_hessian(
      data[nrow(data), ], theta_hat[, i], hessian[, , i], family, G
    )
    gradient <- cost_gradient(data[nrow(data), ], theta_hat[, i], family)
    momentum_step <- solve(hessian[, , i] + epsilon * diag(1, p), gradient)
    momentum <- momentum_coef * momentum - momentum_step
    theta_hat[, i] <- theta_hat[, i] + momentum
    theta_hat[, i] <- DescTools::Winsorize(
      x = theta_hat[, i], minval = winsorise_minval, maxval = winsorise_maxval
    )
  } else if (family == "gaussian") {
    hessian[, , i] <- hessian[, , i] + new_data_x %o% new_data_x
    # B <- as.vector(cmatrix_inv%*%new_data_x)
    # cmatrix_inv <- cmatrix_inv - B%o%B/(1+sum(new_data_x*B))
    gradient <- c(-(new_data_y - new_data_x %*% theta_hat[, i])) * new_data_x
    momentum <- momentum_coef * momentum - solve(hessian[, , i], gradient)
    theta_hat[, i] <- theta_hat[, i] + momentum
    # the choice of norm affects the speed.
    # Spectral norm is more accurate but slower than F norm.
    nc <- norm(hessian[, , i], type = "F")
    theta_hat[, i] <- sign(theta_hat[, i]) * pmax(abs(theta_hat[, i]) - lambda / nc, 0)
  } else {
    stop("family must be one of 'gaussian', 'binomial', or 'poisson'")
  }

  for (kk in 1 + seq_len(k - 1)) {
    for (j in (tau + 1):nrow(data)) {
      if (family == "binomial") {
        prob <- 1 / (1 + exp(-data[j, -1] %*% theta_hat[, i]))
        new_hessian <- (data[j, -1] %o% data[j, -1]) * c((1 - prob) * prob)
        hessian[, , i] <- hessian[, , i] + new_hessian
        gradient <- c(-(data[j, 1] - prob)) * data[j, -1]
        momentum_step <- solve(hessian[, , i] + epsilon * diag(1, p), gradient)
        momentum <- momentum_coef * momentum - momentum_step
        theta_hat[, i] <- theta_hat[, i] + momentum
      } else if (family == "poisson") {
        prob <- exp(data[j, -1] %*% theta_hat[, i])
        new_hessian <- (data[j, -1] %o% data[j, -1]) * min(c(prob), G)
        hessian[, , i] <- hessian[, , i] + new_hessian
        gradient <- c(-(data[j, 1] - prob)) * data[j, -1]
        momentum_step <- solve(hessian[, , i] + epsilon * diag(1, p), gradient)
        momentum <- momentum_coef * momentum - momentum_step
        theta_hat[, i] <- theta_hat[, i] + momentum
        theta_hat[, i] <- DescTools::Winsorize(
          x = theta_hat[, i],
          minval = winsorise_minval,
          maxval = winsorise_maxval
        )
      } else if (family == "gaussian") {
        hessian[, , i] <- hessian[, , i] + data[j, -1] %o% data[j, -1]
        gradient <- c(-(data[j, 1] - data[j, -1] %*% theta_hat[, i])) * data[j, -1]
        momentum <- momentum_coef * momentum - solve(hessian[, , i], gradient)
        theta_hat[, i] <- theta_hat[, i] + momentum
        # the choice of norm affects the speed.
        # Spectral norm is more accurate but slower than F norm.
        nc <- norm(hessian[, , i], type = "F")
        theta_hat[, i] <- sign(theta_hat[, i]) * pmax(abs(theta_hat[, i]) - lambda / nc, 0)
      } else {
        stop("family must be one of 'gaussian', 'binomial', or 'poisson'")
      }
    }
  }

  theta_sum[, i] <- theta_sum[, i] + theta_hat[, i]
  return(list(theta_hat[, i], theta_sum[, i], hessian[, , i], momentum))
}
