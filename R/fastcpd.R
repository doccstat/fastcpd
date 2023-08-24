#' fastcpd: A package for finding change points in an efficient way
#'
#' The fastcpd package provides a function \code{\link{fastcpd}} to find change
#' points in a data set. The function is based on the paper "Sequential Gradient
#' Descent and Quasi-Newton’s Method for Change-Point Analysis" by
#' Xianyang Zhang and Trisha Dawn.
#'
#' @section CITATION:
#' Zhang, Xianyang, and Trisha Dawn.
#' "Sequential Gradient Descent and Quasi-Newton's Method for Change-Point Analysis."
#' arXiv preprint arXiv:2210.12235 (2022).
#'
#' @docType package
#' @name fastcpd
#' @useDynLib fastcpd, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom methods show
NULL
#> NULL


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
#' @param family Family of the models. Can be "binomial", "poisson", "lasso" or
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
#' @param p Number of parameters to be estimated.
#' @param cost Cost function to be used. If not specified, the default is
#'   the negative log-likelihood for the corresponding family.
#' @param cost_gradient Gradient for custom cost function.
#' @param cost_hessian Hessian for custom cost function.
#' @param cp_only Whether to return only the change points or with the cost
#'   values for each segment.
#'
#' @return A class \code{fastcpd} object.
#' @export
#' @examples
#' # Linear regression
#' library(fastcpd)
#' set.seed(1)
#' p <- 3
#' x <- mvtnorm::rmvnorm(300, rep(0, p), diag(p))
#' theta_0 <- rbind(c(1, 1.2, -1), c(-1, 0, 0.5), c(0.5, -0.3, 0.2))
#' y <- c(
#'   x[1:100, ] %*% theta_0[1, ] + rnorm(100, 0, 1),
#'   x[101:200, ] %*% theta_0[2, ] + rnorm(100, 0, 1),
#'   x[201:300, ] %*% theta_0[3, ] + rnorm(100, 0, 1)
#' )
#' result <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data.frame(y = y, x = x),
#'   family = "gaussian"
#' )
#' plot(result)
#' summary(result)
#'
#' # Logistic regression
#' library(fastcpd)
#' set.seed(1)
#' x <- matrix(rnorm(1500, 0, 1), ncol = 5)
#' theta <- rbind(rnorm(5, 0, 1), rnorm(5, 2, 1))
#' y <- c(
#'   rbinom(125, 1, 1 / (1 + exp(-x[1:125, ] %*% theta[1, ]))),
#'   rbinom(175, 1, 1 / (1 + exp(-x[126:300, ] %*% theta[2, ])))
#' )
#' result <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data.frame(y = y, x = x),
#'   family = "binomial"
#' )
#' summary(result)
#'
#' # Poisson regression
#' library(fastcpd)
#' set.seed(1)
#' p <- 3
#' x <- mvtnorm::rmvnorm(1500, rep(0, p), diag(p))
#' delta <- rnorm(p)
#' theta_0 <- c(1, 1.2, -1)
#' y <- c(
#'   rpois(300, exp(x[1:300, ] %*% theta_0)),
#'   rpois(400, exp(x[301:700, ] %*% (theta_0 + delta))),
#'   rpois(300, exp(x[701:1000, ] %*% theta_0)),
#'   rpois(100, exp(x[1001:1100, ] %*% (theta_0 - delta))),
#'   rpois(200, exp(x[1101:1300, ] %*% theta_0)),
#'   rpois(200, exp(x[1301:1500, ] %*% (theta_0 + delta)))
#' )
#' result <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data.frame(y = y, x = x),
#'   beta = (p + 1) * log(1500) / 2,
#'   k = function(x) 0,
#'   family = "poisson",
#'   epsilon = 1e-5
#' )
#' summary(result)
#'
#' # Penalized linear regression
#' library(fastcpd)
#' set.seed(1)
#' n <- 1500
#' p_true <- 6
#' p <- 50
#' x <- mvtnorm::rmvnorm(1500, rep(0, p), diag(p))
#' theta_0 <- rbind(
#'   runif(p_true, -5, -2),
#'   runif(p_true, -3, 3),
#'   runif(p_true, 2, 5),
#'   runif(p_true, -5, 5)
#' )
#' theta_0 <- cbind(theta_0, matrix(0, ncol = p - p_true, nrow = 4))
#' y <- c(
#'   x[1:300, ] %*% theta_0[1, ] + rnorm(300, 0, 1),
#'   x[301:700, ] %*% theta_0[2, ] + rnorm(400, 0, 1),
#'   x[701:1000, ] %*% theta_0[3, ] + rnorm(300, 0, 1),
#'   x[1001:1500, ] %*% theta_0[4, ] + rnorm(500, 0, 1)
#' )
#' result <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data.frame(y = y, x = x),
#'   family = "lasso"
#' )
#' plot(result)
#' summary(result)
#'
#' # Custom cost function: logistic regression
#' library(fastcpd)
#' set.seed(1)
#' p <- 5
#' x <- matrix(rnorm(375 * p, 0, 1), ncol = p)
#' theta <- rbind(rnorm(p, 0, 1), rnorm(p, 2, 1))
#' y <- c(
#'   rbinom(200, 1, 1 / (1 + exp(-x[1:200, ] %*% theta[1, ]))),
#'   rbinom(175, 1, 1 / (1 + exp(-x[126:300, ] %*% theta[2, ])))
#' )
#' data <- data.frame(y = y, x = x)
#' result_builtin <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data,
#'   family = "binomial"
#' )
#' logistic_loss <- function(data, theta) {
#'   x <- data[, -1]
#'   y <- data[, 1]
#'   u <- x %*% theta
#'   nll <- -y * u + log(1 + exp(u))
#'   nll[u > 10] <- -y[u > 10] * u[u > 10] + u[u > 10]
#'   sum(nll)
#' }
#' logistic_loss_gradient <- function(data, theta) {
#'   x <- data[nrow(data), -1]
#'   y <- data[nrow(data), 1]
#'   c(-(y - 1 / (1 + exp(-x %*% theta)))) * x
#' }
#' logistic_loss_hessian <- function(data, theta) {
#'   x <- data[nrow(data), -1]
#'   prob <- 1 / (1 + exp(-x %*% theta))
#'   (x %o% x) * c((1 - prob) * prob)
#' }
#' result_custom <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data,
#'   epsilon = 1e-5,
#'   cost = logistic_loss,
#'   cost_gradient = logistic_loss_gradient,
#'   cost_hessian = logistic_loss_hessian
#' )
#' cat(
#'   "Change points detected by built-in logistic regression model:",
#'   result_builtin@cp_set, "\n",
#'   "Change points detected by custom logistic regression model:",
#'   result_custom@cp_set, "\n",
#'   sep = ""
#' )
#'
#' # Custom cost function: mean shift
#' library(fastcpd)
#' set.seed(1)
#' p <- 1
#' data <- rbind(
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(100, p)),
#'   mvtnorm::rmvnorm(400, mean = rep(50, p), sigma = diag(100, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(2, p), sigma = diag(100, p))
#' )
#' segment_count_guess <- 10
#' block_size <- max(floor(sqrt(nrow(data)) / (segment_count_guess + 1)), 2)
#' block_count <- ceiling(nrow(data) / block_size)
#' data_all_vars <- rep(0, block_count)
#' for (block_index in seq_len(block_count)) {
#'   block_start <- (block_index - 1) * block_size + 1
#'   block_end <- min(block_index * block_size, nrow(data))
#'   data_all_vars[block_index] <- var(data[block_start:block_end, ])
#' }
#' data_all_var <- mean(data_all_vars)
#' mean_loss <- function(data) {
#'   n <- nrow(data)
#'   (norm(data, type = "F")^2 - colSums(data)^2 / n) / 2 / data_all_var +
#'     n / 2 * (log(data_all_var) + log(2 * pi))
#' }
#' mean_loss_result <- fastcpd(
#'   formula = ~ . - 1,
#'   data = data.frame(data),
#'   beta = (p + 1) * log(nrow(data)) / 2,
#'   p = p,
#'   cost = mean_loss
#' )
#' summary(mean_loss_result)
#'
#' # Custom cost function: variance change
#' library(fastcpd)
#' set.seed(1)
#' p <- 1
#' data <- rbind.data.frame(
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(400, mean = rep(0, p), sigma = diag(50, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(2, p))
#' )
#' data_all_mu <- colMeans(data)
#' var_loss <- function(data) {
#'   demeaned_data_norm <- norm(sweep(data, 2, data_all_mu), type = "F")
#'   nrow(data) * (1 + log(2 * pi) + log(demeaned_data_norm^2 / nrow(data))) / 2
#' }
#' var_loss_result <- fastcpd(
#'   formula = ~ . - 1,
#'   data = data,
#'   beta = (p + 1) * log(nrow(data)) / 2,
#'   p = p,
#'   cost = var_loss
#' )
#' summary(var_loss_result)
#'
#' # Custom cost function: mean shift and variance change
#' library(fastcpd)
#' set.seed(1)
#' p <- 1
#' data <- rbind.data.frame(
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(50, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(0, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(400, mean = rep(10, p), sigma = diag(1, p)),
#'   mvtnorm::rmvnorm(300, mean = rep(10, p), sigma = diag(50, p))
#' )
#' meanvar_loss <- function(data) {
#'   loss_part <- (colSums(data^2) - colSums(data)^2 / nrow(data)) / nrow(data)
#'   nrow(data) * (1 + log(2 * pi) + log(loss_part)) / 2
#' }
#' meanvar_loss_result <- fastcpd(
#'   formula = ~ . - 1,
#'   data = data,
#'   beta = (2 * p + 1) * log(nrow(data)) / 2,
#'   p = 2 * p,
#'   cost = meanvar_loss
#' )
#' summary(meanvar_loss_result)
#'
#' # Custom cost function: Huber loss
#' library(fastcpd)
#' set.seed(1)
#' n <- 400 + 300 + 500
#' p <- 5
#' x <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = diag(p))
#' theta <- rbind(
#'   mvtnorm::rmvnorm(1, mean = rep(0, p - 3), sigma = diag(p - 3)),
#'   mvtnorm::rmvnorm(1, mean = rep(5, p - 3), sigma = diag(p - 3)),
#'   mvtnorm::rmvnorm(1, mean = rep(9, p - 3), sigma = diag(p - 3))
#' )
#' theta <- cbind(theta, matrix(0, 3, 3))
#' theta <- theta[rep(seq_len(3), c(400, 300, 500)), ]
#' y_true <- rowSums(x * theta)
#' factor <- c(
#'   2 * stats::rbinom(400, size = 1, prob = 0.95) - 1,
#'   2 * stats::rbinom(300, size = 1, prob = 0.95) - 1,
#'   2 * stats::rbinom(500, size = 1, prob = 0.95) - 1
#' )
#' y <- factor * y_true + stats::rnorm(n)
#' data <- cbind.data.frame(y, x)
#' huber_threshold <- 1
#' huber_loss <- function(data, theta) {
#'   residual <- data[, 1] - data[, -1, drop = FALSE] %*% theta
#'   indicator <- abs(residual) <= huber_threshold
#'   sum(
#'     residual^2 / 2 * indicator +
#'       huber_threshold * (abs(residual) - huber_threshold / 2) * (1 - indicator)
#'   )
#' }
#' huber_loss_gradient <- function(data, theta) {
#'   residual <- c(data[nrow(data), 1] - data[nrow(data), -1] %*% theta)
#'   if (abs(residual) <= huber_threshold) {
#'     -residual * data[nrow(data), -1]
#'   } else {
#'     -huber_threshold * sign(residual) * data[nrow(data), -1]
#'   }
#' }
#' huber_loss_hessian <- function(data, theta) {
#'   residual <- c(data[nrow(data), 1] - data[nrow(data), -1] %*% theta)
#'   if (abs(residual) <= huber_threshold) {
#'     outer(data[nrow(data), -1], data[nrow(data), -1])
#'   } else {
#'     0.01 * diag(length(theta))
#'   }
#' }
#' huber_regression_result <- fastcpd(
#'   formula = y ~ . - 1,
#'   data = data,
#'   beta = (p + 1) * log(n) / 2,
#'   cost = huber_loss,
#'   cost_gradient = huber_loss_gradient,
#'   cost_hessian = huber_loss_hessian
#' )
#' summary(huber_regression_result)
fastcpd <- function(
    formula = y ~ . - 1,
    data,
    beta = NULL,
    segment_count = 10,
    trim = 0.025,
    momentum_coef = 0,
    k = function(x) 0,
    family = NULL,
    epsilon = 1e-10,
    min_prob = 10^10,
    winsorise_minval = -20,
    winsorise_maxval = 20,
    p = NULL,
    cost = negative_log_likelihood,
    cost_gradient = cost_update_gradient,
    cost_hessian = cost_update_hessian,
    cp_only = FALSE) {
  # The following code is adapted from the `lm` function from base R.
  match_formula <- match.call(expand.dots = FALSE)
  matched_formula <- match(c("formula", "data"), names(match_formula), 0L)
  match_formula <- match_formula[c(1L, matched_formula)]
  match_formula$drop.unused.levels <- TRUE
  match_formula[[1L]] <- quote(stats::model.frame)
  match_formula <- eval(match_formula, parent.frame())
  mt <- attr(match_formula, "terms")
  y <- stats::model.response(match_formula, "numeric")
  x <- stats::model.matrix(mt, match_formula)
  data <- cbind(y, x)

  if (is.null(family)) {
    family <- "custom"
  }
  n <- nrow(data)
  if (is.null(p)) {
    p <- ncol(data) - 1
  }
  if (is.null(beta)) {
    beta <- (p + 1) * log(nrow(data)) / 2
  }

  # User provided cost function with explicit expression.
  result <- if (length(formals(cost)) == 1) {
    fastcpd_vanilla_custom(
      data, n, beta, segment_count, trim, momentum_coef, k, epsilon,
      min_prob, winsorise_minval, winsorise_maxval, p, cost, cp_only
    )
    # } else if (vanilla) {
    #   fastcpd_vanilla_custom(
    #     data, n, beta, segment_count, trim, momentum_coef, k, epsilon,
    #     min_prob, winsorise_minval, winsorise_maxval, p,
    #     function(data) {
    #       cost(data = data, theta = NULL, family = family, lambda = 0)
    #     }, cp_only, warm_start
    #   )
  } else {
    fastcpd_builtin(
      data, n, beta, segment_count, trim, momentum_coef, k, family, epsilon,
      min_prob, winsorise_minval, winsorise_maxval, p, cost, cost_gradient,
      cost_hessian, cp_only
    )
  }
  methods::new(
    Class = "fastcpd",
    call = match.call(),
    data = data.frame(data),
    family = family,
    cp_set = result$cp_set,
    cost_values = result$cost_values,
    residuals = result$residual,
    thetas = result$thetas,
    cp_only = cp_only
  )
}

fastcpd_vanilla_custom <- function(
    data, n, beta, segment_count, trim, momentum_coef, k, epsilon,
    min_prob, winsorise_minval, winsorise_maxval, p, cost, cp_only) {
  # fastcpd_vanilla(
  #   data, beta, segment_count, trim, momentum_coef, k, family, epsilon,
  #   min_prob, winsorise_minval, winsorise_maxval, p, cost, cp_only
  # )

  # After t = 1, the r_t_set R_t contains 0 and 1.
  r_t_set <- c(0, 1)
  # C(0)=NULL, C(1)={0}
  cp_set <- append(list(NULL), rep(list(0), n))
  # Objective function: F(0) = -beta
  f_t <- c(-beta, rep(0, n))

  start <- matrix(0, p, n)

  for (t in 2:n) {
    r_t_count <- length(r_t_set)

    # number of cost values is the same as number of elemnts in R_t
    cval <- rep(0, r_t_count)

    # for tau in R_t\{t-1}
    for (i in 1:(r_t_count - 1)) {
      tau <- r_t_set[i]

      if (t - tau >= 1) {
        # if (warm_start && t - tau >= 50) {
        #   cost_result <- cost(data[(tau + 1):t, , drop = FALSE], start = start[, tau + 1])
        #   start[, tau + 1] <- cost_result$par
        #   cval[i] <- cost_result$value
        # } else {
        cval[i] <- cost(data[(tau + 1):t, , drop = FALSE])
        # }
      }
    }

    # Step 3
    cval[r_t_count] <- 0
    obj <- cval + f_t[r_t_set + 1] + beta
    min_val <- min(obj)
    tau_star <- r_t_set[which(obj == min_val)[1]]

    # Step 4
    cp_set[[t + 1]] <- c(cp_set[[tau_star + 1]], tau_star)
    # print(r_t_set)
    # print(cp_set[[t + 1]])
    # print(cval)
    # print(f_t[r_t_set + 1])
    # print(obj)
    # print(min_val)
    # print((cval + f_t[r_t_set + 1]) <= min_val)

    # Step 5
    pruned_left <- (cval + f_t[r_t_set + 1]) <= min_val
    r_t_set <- c(r_t_set[pruned_left], t)

    # Objective function F(t).
    f_t[t + 1] <- min_val
  }
  # print(cp_set)

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
  thetas <- matrix(NA, nrow = 0, ncol = 0)

  if (cp_only) {
    cost_values <- numeric(0)
  } else {
    cp_loc <- unique(c(0, cp_set, n))
    cost_values <- rep(0, length(cp_loc) - 1)
    for (i in 1:(length(cp_loc) - 1)) {
      cost_values[i] <- cost(data[(cp_loc[i] + 1):cp_loc[i + 1], , drop = FALSE])
    }
  }
  thetas <- data.frame(thetas)
  list(
    cp_set = cp_set,
    cost_values = cost_values,
    residual = 0,
    thetas = thetas
  )
}

fastcpd_builtin <- function(
    data, n, beta, segment_count, trim, momentum_coef, k, family, epsilon,
    min_prob, winsorise_minval, winsorise_maxval, p, cost, cost_gradient,
    cost_hessian, cp_only) {
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
    segment_theta <- if (family == "custom") {
      if (p == 1) {
        optim_result <- stats::optim(
          par = 0,
          fn = function(theta, data) {
            cost(data = data, theta = log(theta / (1 - theta)))
          },
          method = "Brent",
          lower = 0,
          upper = 1,
          data = data_segment
        )
        log(optim_result$par / (1 - optim_result$par))
      } else {
        stats::optim(
          par = rep(0, p),
          fn = cost,
          data = data_segment,
          method = "L-BFGS-B"
        )$par
      }
    } else {
      cost(
        data = data_segment,
        theta = NULL,
        family = family,
        lambda = 0,
        cv = TRUE
      )$par
    }
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
    theta_sum <- theta_hat <- DescTools::Winsorize(
      matrix(segment_theta_hat[1, ]),
      minval = winsorise_minval,
      maxval = winsorise_maxval
    )
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
  } else if (family == "custom") {
    theta_sum <- theta_hat <- matrix(segment_theta_hat[1, ])
    hessian <- array(
      matrix(0, p, p),
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

      cost_update_result <- cost_update(
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
        lambda = lambda,
        cost_gradient = cost_gradient,
        cost_hessian = cost_hessian
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
        cval[i] <- cost(data[(tau + 1):t, , drop = FALSE], theta, family, lambda)$value
      } else if (family == "custom" && t - tau >= 1) {
        cval[i] <- cost(data[(tau + 1):t, , drop = FALSE], theta)
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
    } else if (family == "custom") {
      cum_coef_add <- coef_add <- segment_theta_hat[segment_indices[t], ]
      hessian_new <- matrix(0, p, p)
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
  residual <- 0[family == "custom"]

  if (cp_only) {
    cost_values <- numeric(0)
    thetas <- matrix(NA, nrow = 0, ncol = 0)
  } else {
    cp_loc <- unique(c(0, cp_set, n))
    for (i in 1:(length(cp_loc) - 1)) {
      if (family == "custom") {
        if (p == 1) {
          optim_result <- stats::optim(
            par = 0,
            fn = function(theta, data) {
              cost(data = data, theta = log(theta / (1 - theta)))
            },
            method = "Brent",
            lower = 0,
            upper = 1,
            data = data[(cp_loc[i] + 1):cp_loc[i + 1], , drop = FALSE]
          )
          cost_result <- list(
            par = log(optim_result$par / (1 - optim_result$par)),
            value = exp(optim_result$value) / (1 + exp(optim_result$value))
          )
        } else {
          cost_result <- stats::optim(
            par = rep(0, p),
            fn = cost,
            data = data[(cp_loc[i] + 1):cp_loc[i + 1], , drop = FALSE],
            method = "L-BFGS-B"
          )
        }
      } else {
        cost_result <- cost(data[(cp_loc[i] + 1):cp_loc[i + 1], , drop = FALSE], NULL, family, lambda)
        residual <- c(residual, cost_result$residuals)
      }
      cost_values <- c(cost_values, cost_result$value)
      thetas <- cbind(thetas, cost_result$par)
    }
  }
  thetas <- data.frame(thetas)
  if (ncol(thetas) > 0) {
    names(thetas) <- paste0("segment ", seq_len(ncol(thetas)))
  }
  list(
    cp_set = cp_set,
    cost_values = cost_values,
    residual = residual,
    thetas = thetas
  )
}
