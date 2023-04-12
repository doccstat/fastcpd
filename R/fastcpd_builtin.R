#' Sequential Gradient Descent and Quasi-Newton’s Method for Change-Point
#' Analysis
#'
#' @param formula A symbolic description of the model to be fitted.
#' @param data A data frame containing the data to be segmented.
#' @param beta Initial cost value.
#' @param segment_count Number of segments for initial guess.
#' @param trim Trimming for the boundary change points.
#' @param momentum_coef Momentum coefficient to be applied to each update.
#' @param k Function on number of epochs in SGD.
#' @param family Family of the model. Can be "binomial", "poisson", "lasso" or
#'   "gaussian". If not provided, the user must specify the cost function and
#'   its gradient (and Hessian).
#' @param epsilon Epsilon to avoid numerical issues. Only used for binomial and
#'   poisson.
#' @param min_prob Minimum probability to avoid numerical issues. Only used for
#'   poisson.
#' @param winsorise_minval Minimum value to be winsorised. Only used for
#'   poisson.
#' @param winsorise_maxval Maximum value to be winsorised. Only used for
#'   poisson.
#' @param cp_only Whether to return only the change points or with the cost
#'   values for each segment.
#'
#' @return Change points and corresponding cost values.
fastcpd_builtin <- function(
  formula, data, beta, segment_count, trim, momentum_coef, k, family,
  epsilon, min_prob, winsorise_minval, winsorise_maxval, p, cp_only
) {
  n <- nrow(data)
  if (is.null(p)) {
    p <- ncol(data) - 1
  }
  if (is.null(beta)) {
    beta <- (p + 1) * log(nrow(data)) / 2
  }

  if (family %in% c("lasso", "gaussian")) {
    err_sd <- act_num <- rep(NA, segment_count)
  }

  # After t = 1, the r_t_set R_t contains 0 and 1.
  r_t_set <- c(0, 1)
  # C(0)=NULL, C(1)={0}
  cp_set <- append(list(NULL), rep(list(0), n))
  # Objective function: F(0) = -beta
  f_t <- c(-beta, rep(0, n))
  momentum <- rep(0, p)

  # choose the initial values based on pre-segmentation

  segment_indices <- ceiling(seq_len(n) / ceiling(n / segment_count))
  segment_theta_hat <- matrix(NA, segment_count, p)
  # Remark 3.4: initialize theta_hat_t_t to be the estimate in the segment
  for (segment_index in seq_len(segment_count)) {
    data_segment <- data[segment_indices == segment_index, , drop = FALSE]
    segment_theta <- negative_log_likelihood(
      data = data_segment,
      theta = NULL,
      family = family,
      lambda = NULL,
      cv = TRUE
    )$par
    segment_theta_hat[segment_index, ] <- segment_theta
    if (family %in% c("lasso", "gaussian")) {
      response_estimate <- data_segment[, -1, drop = FALSE] %*% c(segment_theta)
      segment_residual <- data_segment[, 1] - response_estimate
      err_sd[segment_index] <- sqrt(mean(segment_residual^2))
      act_num[segment_index] <- sum(abs(segment_theta) > 0)
    }
  }

  if (family %in% c("lasso", "gaussian")) {
    # only works if error sd is unchanged.
    err_sd_mean <- mean(err_sd)

    act_num_mean <- mean(act_num)

    # seems to work but there might be better choices
    beta <- (act_num_mean + 1) * beta
  }

  # For the first data.
  if (family == "binomial") {
    theta_sum <- theta_hat <- matrix(segment_theta_hat[1, ])
    prob <- 1 / (1 + exp(-theta_hat %*% data[1, -1]))
    hessian <- array(
      (data[1, -1] %o% data[1, -1]) * c(prob * (1 - prob)),
      c(p, p, 1)
    )
  } else if (family == "poisson") {
    theta_sum <- theta_hat <- matrix(DescTools::Winsorize(
      matrix(segment_theta_hat[1, ]),
      minval = winsorise_minval,
      maxval = winsorise_maxval
    ))
    hessian <- array(
      (data[1, -1] %o% data[1, -1]) * c(exp(theta_hat %*% data[1, -1])),
      c(p, p, 1)
    )
  } else if (family %in% c("lasso", "gaussian")) {
    theta_sum <- theta_hat <- matrix(segment_theta_hat[1, ])
    hessian <- array(
      data[1, -1] %o% data[1, -1] + epsilon * diag(1, p),
      c(p, p, 1)
    )
  }

  for (t in 2:n) {
    r_t_count <- length(r_t_set)
    # number of cost values is the same as number of elemnts in R_t
    cval <- rep(0, r_t_count)

    # for tau in R_t\{t-1}
    for (i in 1:(r_t_count - 1)) {
      tau <- r_t_set[i]
      if (family == "lasso") {
        lambda <- err_sd_mean * sqrt(2 * log(p) / (t - tau))
      } else {
        lambda <- 0
      }

      cost_update_result <- cost_update_builtin(
        data = data[seq_len(t), , drop = FALSE],
        theta_hat = theta_hat,
        theta_sum = theta_sum,
        hessian = hessian,
        tau = tau,
        i = i,
        k = k,
        family = family,
        momentum = momentum,
        momentum_coef = momentum_coef,
        min_prob = min_prob,
        winsorise_minval = winsorise_minval,
        winsorise_maxval = winsorise_maxval,
        epsilon = epsilon,
        lambda = lambda
      )
      theta_hat[, i] <- cost_update_result[[1]]
      theta_sum[, i] <- cost_update_result[[2]]
      hessian[, , i] <- cost_update_result[[3]]
      momentum <- cost_update_result[[4]]

      tau <- r_t_set[i]
      theta <- theta_sum[, i] / (t - tau)
      if (family == "poisson" && t - tau >= p) {
        theta <- DescTools::Winsorize(
          theta,
          minval = winsorise_minval,
          maxval = winsorise_maxval
        )
      }

      if (
        (family %in% c("binomial", "poisson") && t - tau >= p) ||
        (family %in% c("lasso", "gaussian") && t - tau >= 3)
      ) {
        cval[i] <- negative_log_likelihood(data[(tau + 1):t, , drop = FALSE], theta, family, lambda)
      }
    }

    # the choice of initial values requires further investigation

    # for tau = t-1
    new_data <- data[t, -1]
    if (family == "binomial") {
      cum_coef_add <- coef_add <- segment_theta_hat[segment_indices[t], ]
      prob <- 1 / (1 + exp(-coef_add %*% new_data))
      hessian_new <- (new_data %o% new_data) * c(prob * (1 - prob))
    } else if (family == "poisson") {
      cum_coef_add <- coef_add <- DescTools::Winsorize(
        x = segment_theta_hat[segment_indices[t], ],
        minval = winsorise_minval,
        maxval = winsorise_maxval
      )
      hessian_new <- (new_data %o% new_data) * c(exp(coef_add %*% new_data))
    } else if (family %in% c("lasso", "gaussian")) {
      cum_coef_add <- coef_add <- segment_theta_hat[segment_indices[t], ]
      hessian_new <- new_data %o% new_data + epsilon * diag(1, p)
    }

    theta_hat <- cbind(theta_hat, coef_add)
    theta_sum <- cbind(theta_sum, cum_coef_add)
    hessian <- abind::abind(hessian, hessian_new, along = 3)

    # Step 3
    cval[r_t_count] <- 0
    obj <- cval + f_t[r_t_set + 1] + beta
    min_val <- min(obj)
    tau_star <- r_t_set[which(obj == min_val)[1]]

    # Step 4
    cp_set[[t + 1]] <- c(cp_set[[tau_star + 1]], tau_star)

    # Step 5
    pruned_left <- (cval + f_t[r_t_set + 1]) <= min_val
    r_t_set <- c(r_t_set[pruned_left], t)

    theta_hat <- theta_hat[, pruned_left, drop = FALSE]
    theta_sum <- theta_sum[, pruned_left, drop = FALSE]
    hessian <- hessian[, , pruned_left, drop = FALSE]

    # Objective function F(t).
    f_t[t + 1] <- min_val
  }

  # Remove change-points close to the boundaries

  cp_set <- cp_set[[n + 1]]
  cp_set <- cp_set[(cp_set >= trim * n) & (cp_set <= (1 - trim) * n)]
  cp_set <- sort(unique(c(0, cp_set)))

  segment_indices <- which((diff(cp_set) < trim * n) == TRUE)
  if (length(segment_indices) > 0) {
    cp_set <- floor(
      (cp_set[-(segment_indices + 1)] + cp_set[-segment_indices]) / 2
    )
  }
  cp_set <- cp_set[cp_set > 0]
  cost_values <- thetas <- NULL
  residual <- NULL

  if (cp_only) {
    cost_values <- numeric(0)
    thetas <- matrix(NA, nrow = 0, ncol = 0)
  } else {
    cp_loc <- unique(c(0, cp_set, n))
    for (i in 1:(length(cp_loc) - 1)) {
      cost_result <- negative_log_likelihood(data[(cp_loc[i] + 1):cp_loc[i + 1], , drop = FALSE], NULL, family, lambda)
      residual <- c(residual, cost_result$residuals)
      cost_values <- c(cost_values, cost_result$value)
      thetas <- cbind(thetas, cost_result$par)
    }
  }

  methods::new(
    Class = "fastcpd",
    call = match.call(),
    data = data.frame(data),
    family = family,
    cp_set = cp_set,
    cost_values = cost_values,
    residuals = residual,
    thetas = thetas,
    cp_only = cp_only
  )
}
