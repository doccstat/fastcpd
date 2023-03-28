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
#' @param epsilon Epsilon to avoid numerical issues. Only used for binomial and
#'   poisson.
#' @param min_prob Minimum probability to avoid numerical issues. Only used for
#'   poisson.
#' @param winsorise_minval Minimum value to be winsorised. Only used for
#'   poisson.
#' @param winsorise_maxval Maximum value to be winsorised. Only used for
#'   poisson.
#' @param lambda Lambda for L1 regularization. Only used for lasso.
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
  epsilon,
  min_prob,
  winsorise_minval,
  winsorise_maxval,
  lambda,
  cost_gradient,
  cost_hessian
) {
  hessian[, , i] <- if (family == "custom") {
    cost_hessian(data[nrow(data), ], theta_hat[, i], hessian[, , i])
  } else {
    cost_hessian(
      data[nrow(data), ], theta_hat[, i], hessian[, , i], family, min_prob
    )
  }
  gradient <- if (family == "custom") {
    cost_gradient(data[nrow(data), ], theta_hat[, i])
  } else {
    cost_gradient(data[nrow(data), ], theta_hat[, i], family)
  }
  hessian_psd <- hessian[, , i] + epsilon * diag(1, nrow(theta_hat))
  momentum_step <- solve(hessian_psd, gradient)
  momentum <- momentum_coef * momentum - momentum_step
  theta_hat[, i] <- theta_hat[, i] + momentum

  if (family == "poisson") {
    theta_hat[, i] <- DescTools::Winsorize(
      x = theta_hat[, i], minval = winsorise_minval, maxval = winsorise_maxval
    )
  } else if (family %in% c("lasso", "gaussian")) {
    # the choice of norm affects the speed.
    # Spectral norm is more accurate but slower than F norm.
    hessian_norm <- norm(as.matrix(hessian[, , i]), type = "F")
    normd <- abs(theta_hat[, i]) - lambda / hessian_norm
    theta_hat[, i] <- sign(theta_hat[, i]) * pmax(normd, 0)
  }

  for (kk in 1 + seq_len(k(nrow(data) - tau))) {
    for (j in (tau + 1):nrow(data)) {
      hessian[, , i] <- if (family == "custom") {
        cost_hessian(data[j, ], theta_hat[, i], hessian[, , i])
      } else {
        cost_hessian(
          data[j, ], theta_hat[, i], hessian[, , i], family, min_prob
        )
      }
      gradient <- if (family == "custom") {
        cost_gradient(data[j, ], theta_hat[, i])
      } else {
        cost_gradient(data[j, ], theta_hat[, i], family)
      }
      hessian_psd <- hessian[, , i] + epsilon * diag(1, nrow(theta_hat))
      momentum_step <- solve(hessian_psd, gradient)
      momentum <- momentum_coef * momentum - momentum_step
      theta_hat[, i] <- theta_hat[, i] + momentum

      if (family == "poisson") {
        theta_hat[, i] <- DescTools::Winsorize(
          x = theta_hat[, i],
          minval = winsorise_minval,
          maxval = winsorise_maxval
        )
      } else if (family %in% c("lasso", "gaussian")) {
        hessian_norm <- norm(as.matrix(hessian[, , i]), type = "F")
        normd <- abs(theta_hat[, i]) - lambda / hessian_norm
        theta_hat[, i] <- sign(theta_hat[, i]) * pmax(normd, 0)
      }
    }
  }

  theta_sum[, i] <- theta_sum[, i] + theta_hat[, i]
  return(list(theta_hat[, i], theta_sum[, i], hessian[, , i], momentum))
}
