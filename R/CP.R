#' Dynamic Programming with Pruning and Gradient Descent
#'
#' @param data TODO
#' @param beta TODO
#' @param segment_count TODO
#' @param trim TODO
#' @param family TODO
#' @param ... TODO
#'
#' @return TODO
#' @export
CP <- function(data, beta, segment_count = 10, trim = 0.025, family, ...) {
  args_list <- list(...)
  n <- dim(data)[1]
  p <- dim(data)[2] - 1

  # choose the initial values based on pre-segmentation

  index <- rep(1:segment_count, rep(n / segment_count, segment_count))
  theta <- matrix(NA, segment_count, p)
  for (segment_index in 1:segment_count) {
    if (family %in% c("binomial", "poisson")) {
      theta[segment_index, ] <- fastglm::fastglm(
        data[index == segment_index, 1:p, drop = FALSE],
        data[index == segment_index, p + 1],
        family
      )$coefficients
    } else if (family == "gaussian") {
      cvfit <- glmnet::cv.glmnet(as.matrix(data[index == segment_index, 1:p]), data[index == segment_index, p + 1], family = family)
      theta[segment_index, ] <- theta_hat(cvfit, s = "lambda.1se")[-1]
      resi <- data[index == segment_index, p + 1] - as.matrix(data[index == segment_index, 1:p]) %*% as.numeric(theta[segment_index, ])
      err_sd[segment_index] <- sqrt(mean(resi^2))
      act_num[segment_index] <- sum(abs(theta[segment_index, ]) > 0)
    }
  }

  if (family == "gaussian") {
    err_sd_mean <- mean(err_sd) # only works if error sd is unchanged.
    act_num_mean <- mean(act_num)
    beta <- (act_num_mean + 1) * beta # seems to work but there might be better choices
  }

  # t = 1
  if (family == "binomial") {
    theta_sum <- theta_hat <- matrix(theta[1, ], p, 1)
    p <- 1 / (1 + exp(-theta_hat %*% data[1, 1:p]))
    hessian <- array((data[1, 1:p] %o% data[1, 1:p]) * as.numeric(p * (1 - p)), c(p, p, 1))
  } else if (family == "poisson") {
    theta_sum <- theta_hat <- DescTools::Winsorize(matrix(theta[1, ], p, 1), minval = args_list$L, maxval = args_list$H)
    e_eta <- exp(theta_hat %*% data[1, 1:p])
    const <- e_eta
    hessian <- array((data[1, 1:p] %o% data[1, 1:p]) * as.numeric(const), c(p, p, 1))
  } else if (family == "gaussian") {
    theta_sum <- theta_hat <- matrix(theta[1, ], p, 1)
    # eta <- theta_hat %*% data[1, 1:p]
    # c_int <- diag(1/epsilon,p) - data[1, 1:p]%o%data[1, 1:p]/epsilon^2/(1+sum(data[1, 1:p]^2)/epsilon)
    # cmatrix_inv <- array(c_int, c(p,p,1))
    hessian <- array(data[1, 1:p] %o% data[1, 1:p] + args_list$epsilon * diag(1, p), c(p, p, 1))
  }

  # After t = 1, the r_t_set R_t contains 0 and 1.
  r_t_set <- c(0, 1)
  # C(0)=NULL, C(1)={0}
  cp_set <- append(list(NULL), rep(list(0), n))
  # F(0)=-beta
  f_t <- c(-beta, rep(0, n))

  for (t in 2:n) {
    m <- length(r_t_set)
    # number of cost values is the same as number of elemnts in R_t
    cval <- rep(0, m)

    # for tau in R_t\{t-1}
    for (i in 1:(m - 1)) {
      tau <- r_t_set[i]
      if (family == "binomial") {
        out <- cost_logistic_update(data[t, ], theta_hat[, i], theta_sum[, i], hessian[, , i])
        theta_hat[, i] <- out[[1]]
        theta_sum[, i] <- out[[2]]
        hessian[, , i] <- out[[3]]
        if (t - tau >= p) {
          cval[i] <- neg_log_lik(data[(tau + 1):t, ], theta_sum[, i] / (t - tau), family = family)
        }
      } else if (family == "poisson") {
        out <- cost_poisson_update(data[t, ], theta_hat[, i], theta_sum[, i], hessian[, , i], epsilon = args_list$epsilon, G = args_list$G, L = args_list$L, H = args_list$H)
        theta_hat[, i] <- out[[1]]
        theta_sum[, i] <- out[[2]]
        hessian[, , i] <- out[[3]]
        cum_coef_win <- DescTools::Winsorize(theta_sum[, i] / (t - tau), minval = args_list$L, maxval = args_list$H)
        if (t - tau >= p) {
          cval[i] <- neg_log_lik(data[(tau + 1):t, ], cum_coef_win, family = family)
        }
      } else if (family == "gaussian") {
        out <- cost_lasso_update(data[t, ], theta_hat[, i], theta_sum[, i], hessian[, , i], lambda = err_sd_mean * sqrt(2 * log(p) / (t - tau)))
        theta_hat[, i] <- out[[1]]
        theta_sum[, i] <- out[[2]]
        # cmatrix_inv[,,i] <- out[[3]]
        hessian[, , i] <- out[[3]]
        if (t - tau >= 3) {
          cval[i] <- neg_log_lik(data[(tau + 1):t, ], theta_sum[, i] / (t - tau), family = "gaussian", lambda = err_sd_mean * sqrt(2 * log(p) / (t - tau)))
        }
      }
    }

    # the choice of initial values requires further investigation

    new_data <- data[t, 1:p]
    if (family == "binomial") {
      cum_coef_add <- coef_add <- theta[index[t], ]
      p <- 1 / (1 + exp(-coef_add %*% new_data))
      hessian_new <- (new_data %o% new_data) * as.numeric(p * (1 - p))
    } else if (family == "poisson") {
      cum_coef_add <- coef_add <- DescTools::Winsorize(theta[index[t], ], minval = args_list$L, maxval = args_list$H) ####
      hessian_new <- (new_data %o% new_data) * as.numeric(exp(coef_add %*% new_data))
    } else if (family == "gaussian") {
      cum_coef_add <- coef_add <- theta[index[t], ]
      hessian_new <- new_data %o% new_data + args_list$epsilon * diag(1, p)
    }

    theta_hat <- cbind(theta_hat, coef_add)
    theta_sum <- cbind(theta_sum, cum_coef_add)
    hessian <- abind::abind(hessian, hessian_new, along = 3)

    # Adding a momentum term (TBD)

    obj <- cval + f_t[r_t_set + 1] + beta
    min_val <- min(obj)
    ind <- which(obj == min_val)[1]
    tau_star <- r_t_set[ind]
    cp_set[[t + 1]] <- c(cp_set[[tau_star + 1]], tau_star)

    # Step 5
    ind2 <- (cval + f_t[r_t_set + 1]) <= min_val
    r_t_set <- c(r_t_set[ind2], t)

    theta_hat <- theta_hat[, ind2, drop = FALSE]
    theta_sum <- theta_sum[, ind2, drop = FALSE]
    hessian <- hessian[, , ind2, drop = FALSE]
    # F(t)
    f_t[t + 1] <- min_val
  }

  # Remove change-points close to the boundaries

  cp <- cp_set[[n + 1]]
  if (family == "binomial") {
    if (length(cp) > 0) {
      ind3 <- (seq_len(length(cp)))[(cp < trim * n) | (cp > (1 - trim) * n)]
      cp <- cp[-ind3]
    }

    nLL <- 0
    cp_loc <- unique(c(0, cp, n))
    for (i in 1:(length(cp_loc) - 1))
    {
      seg <- (cp_loc[i] + 1):cp_loc[i + 1]
      data_seg <- data[seg, ]
      out <- fastglm::fastglm(
        as.matrix(data_seg[, 1:p]),
        data_seg[, p + 1],
        family
      )
      nLL <- out$deviance / 2 + nLL
    }

    output <- list(cp, nLL)
    names(output) <- c("cp", "nLL")
  } else if (family == "poisson") {
    if (length(cp) > 0) {
      ind3 <- seq_len(length(cp))[(cp < trim * n) | (cp > (1 - trim) * n)]
      if (length(ind3) > 0) cp <- cp[-ind3]
    }

    cp <- sort(unique(c(0, cp)))
    index <- which((diff(cp) < trim * n) == TRUE)
    if (length(index) > 0) cp <- floor((cp[-(index + 1)] + cp[-index]) / 2)
    cp <- cp[cp > 0]

    # nLL <- 0
    # cp_loc <- unique(c(0,cp,n))
    # for(i in 1:(length(cp_loc)-1))
    # {
    #   seg <- (cp_loc[i]+1):cp_loc[i+1]
    #   data_seg <- data[seg,]
    #   out <- fastglm(as.matrix(data_seg[, 1:p]), data_seg[, p+1], family="Poisson")
    #   nLL <- out$deviance/2 + nLL
    # }

    # output <- list(cp, nLL)
    # names(output) <- c("cp", "nLL")

    output <- list(cp)
    names(output) <- c("cp")
  } else if (family == "gaussian") {
    if (length(cp) > 0) {
      ind3 <- seq_len(length(cp))[(cp < trim * n) | (cp > (1 - trim) * n)]
      if (length(ind3) > 0) cp <- cp[-ind3]
    }

    cp <- sort(unique(c(0, cp)))
    index <- which((diff(cp) < trim * n) == TRUE)
    if (length(index) > 0) cp <- floor((cp[-(index + 1)] + cp[-index]) / 2)
    cp <- cp[cp > 0]

    # nLL <- 0
    # cp_loc <- unique(c(0,cp,n))
    # for(i in 1:(length(cp_loc)-1))
    # {
    #  seg <- (cp_loc[i]+1):cp_loc[i+1]
    #  data_seg <- data[seg,]
    #  out <- fastglm(as.matrix(data_seg[, 1:p]), data_seg[, p+1], family="binomial")
    #  nLL <- out$deviance/2 + nLL
    # }

    output <- list(cp)
    names(output) <- c("cp")
  }

  return(output)
}
