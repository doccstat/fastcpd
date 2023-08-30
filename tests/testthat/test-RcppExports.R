run_cpp_tests("fastcpd")

# Everything in this script is provided as is. The purpose of this script is to
# do a sanity check on the C++ implementation of `fastcpd`.

#' Cost function for Logistic regression, i.e. binomial family in GLM.
#'
#' @param data Data to be used to calculate the cost values. The last column is
#'     the response variable.
#' @param family Family of the distribution.
#' @keywords internal
#'
#' @noRd
#' @return Cost value for the corresponding segment of data.
cost_glm_binomial <- function(data, family = "binomial") {
  data <- as.matrix(data)
  p <- dim(data)[2] - 1
  out <- fastglm::fastglm(as.matrix(data[, 1:p]), data[, p + 1], family = family)
  return(out$deviance / 2)
}

#' Implementation of vanilla PELT for logistic regression type data.
#'
#' @param data Data to be used for change point detection.
#' @param beta Penalty coefficient for the number of change points.
#' @param cost Cost function to be used to calculate cost values.
#' @keywords internal
#'
#' @noRd
#' @return A list consisting of two: change point locations and negative log
#'     likelihood values for each segment.
pelt_vanilla_binomial <- function(data, beta, cost = cost_glm_binomial) {
  n <- dim(data)[1]
  p <- dim(data)[2] - 1
  Fobj <- c(-beta, 0)
  cp_set <- list(NULL, 0)
  set <- c(0, 1)
  for (t in 2:n)
  {
    m <- length(set)
    cval <- rep(NA, m)
    for (i in 1:m)
    {
      k <- set[i] + 1
      if (t - k >= p - 1) cval[i] <- suppressWarnings(cost(data[k:t, ])) else cval[i] <- 0
    }
    obj <- cval + Fobj[set + 1] + beta
    min_val <- min(obj)
    ind <- which(obj == min_val)[1]
    cp_set_add <- c(cp_set[[set[ind] + 1]], set[ind])
    cp_set <- append(cp_set, list(cp_set_add))
    ind2 <- (cval + Fobj[set + 1]) <= min_val
    set <- c(set[ind2], t)
    Fobj <- c(Fobj, min_val)
  }
  cp <- cp_set[[n + 1]]
  nLL <- 0
  cp_loc <- unique(c(0, cp, n))
  for (i in 1:(length(cp_loc) - 1))
  {
    seg <- (cp_loc[i] + 1):cp_loc[i + 1]
    data_seg <- data[seg, ]
    out <- fastglm::fastglm(as.matrix(data_seg[, 1:p]), data_seg[, p + 1], family = "binomial")
    nLL <- out$deviance / 2 + nLL
  }

  output <- list(cp, nLL)
  names(output) <- c("cp", "nLL")
  return(output)
}

#' Function to update the coefficients using gradient descent.
#'
#' @param data_new New data point used to update the coeffient.
#' @param coef Previous coeffient to be updated.
#' @param cum_coef Summation of all the past coefficients to be used in
#'     averaging.
#' @param cmatrix Hessian matrix in gradient descent.
#' @param epsilon Small adjustment to avoid singularity when doing inverse on
#'     the Hessian matrix.
#' @keywords internal
#'
#' @noRd
#' @return A list of values containing the new coefficients, summation of
#'     coefficients so far and all the Hessian matrices.
cost_logistic_update <- function(data_new, coef, cum_coef, cmatrix, epsilon = 1e-10) {
  p <- length(data_new) - 1
  X_new <- data_new[1:p]
  Y_new <- data_new[p + 1]
  eta <- X_new %*% coef
  mu <- 1 / (1 + exp(-eta))
  cmatrix <- cmatrix + (X_new %o% X_new) * as.numeric((1 - mu) * mu)
  lik_dev <- as.numeric(-(Y_new - mu)) * X_new
  coef <- coef - solve(cmatrix + epsilon * diag(1, p), lik_dev)
  cum_coef <- cum_coef + coef
  return(list(coef, cum_coef, cmatrix))
}

#' Calculate negative log likelihood given data segment and guess of coefficient.
#'
#' @param data Data to be used to calculate the negative log likelihood.
#' @param b Guess of the coefficient.
#' @keywords internal
#'
#' @noRd
#' @return Negative log likelihood.
neg_log_lik_binomial <- function(data, b) {
  p <- dim(data)[2] - 1
  X <- data[, 1:p, drop = FALSE]
  Y <- data[, p + 1, drop = FALSE]
  u <- as.numeric(X %*% b)
  L <- -Y * u + log(1 + exp(u))
  return(sum(L))
}

#' Find change points using dynamic programming with pruning and SeGD.
#'
#' @param data Data used to find change points.
#' @param beta Penalty coefficient for the number of change points.
#' @param B Initial guess on the number of change points.
#' @param trim Propotion of the data to ignore the change points at the
#'     beginning, ending and between change points.
#' @keywords internal
#'
#' @noRd
#' @return A list containing potential change point locations and negative log
#'     likelihood for each segment based on the change points guess.
segd_binomial <- function(data, beta, B = 10, trim = 0.025) {
  n <- dim(data)[1]
  p <- dim(data)[2] - 1
  Fobj <- c(-beta, 0)
  cp_set <- list(NULL, 0)
  set <- c(0, 1)

  # choose the initial values based on pre-segmentation

  index <- rep(1:B, rep(n / B, B))
  coef.int <- matrix(NA, B, p)
  for (i in 1:B)
  {
    out <- fastglm::fastglm(as.matrix(data[index == i, 1:p]), data[index == i, p + 1], family = "binomial")
    coef.int[i, ] <- coef(out)
  }
  X1 <- data[1, 1:p]
  cum_coef <- coef <- matrix(coef.int[1, ], p, 1)
  e_eta <- exp(coef %*% X1)
  const <- e_eta / (1 + e_eta)^2
  cmatrix <- array((X1 %o% X1) * as.numeric(const), c(p, p, 1))

  for (t in 2:n)
  {
    m <- length(set)
    cval <- rep(NA, m)

    for (i in 1:(m - 1))
    {
      coef_c <- coef[, i]
      cum_coef_c <- cum_coef[, i]
      cmatrix_c <- cmatrix[, , i]
      out <- cost_logistic_update(data[t, ], coef_c, cum_coef_c, cmatrix_c)
      coef[, i] <- out[[1]]
      cum_coef[, i] <- out[[2]]
      cmatrix[, , i] <- out[[3]]
      k <- set[i] + 1
      if (t - k >= p - 1) cval[i] <- neg_log_lik_binomial(data[k:t, ], cum_coef[, i] / (t - k + 1)) else cval[i] <- 0
    }

    # the choice of initial values requires further investigation

    cval[m] <- 0
    Xt <- data[t, 1:p]
    cum_coef_add <- coef_add <- coef.int[index[t], ]
    e_eta_t <- exp(coef_add %*% Xt)
    const <- e_eta_t / (1 + e_eta_t)^2
    cmatrix_add <- (Xt %o% Xt) * as.numeric(const)

    coef <- cbind(coef, coef_add)
    cum_coef <- cbind(cum_coef, cum_coef_add)
    cmatrix <- abind::abind(cmatrix, cmatrix_add, along = 3)

    # Adding a momentum term (TBD)

    obj <- cval + Fobj[set + 1] + beta
    min_val <- min(obj)
    ind <- which(obj == min_val)[1]
    cp_set_add <- c(cp_set[[set[ind] + 1]], set[ind])
    cp_set <- append(cp_set, list(cp_set_add))
    ind2 <- (cval + Fobj[set + 1]) <= min_val
    set <- c(set[ind2], t)
    coef <- coef[, ind2, drop = FALSE]
    cum_coef <- cum_coef[, ind2, drop = FALSE]
    cmatrix <- cmatrix[, , ind2, drop = FALSE]
    Fobj <- c(Fobj, min_val)
  }

  # Remove change-points close to the boundaries

  cp <- cp_set[[n + 1]]
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
    out <- fastglm::fastglm(as.matrix(data_seg[, 1:p]), data_seg[, p + 1], family = "binomial")
    nLL <- out$deviance / 2 + nLL
  }

  output <- list(cp, nLL)
  names(output) <- c("cp", "nLL")
  return(output)
}

test_that("logistic regression", {
  # This is the same example with `fastcpd` documentation. Please keep it in
  # sync if the documentation ever changes.
  set.seed(1)

  kChangePointLocation <- 125
  kNumberOfDataPoints <- 300
  kDimension <- 5

  # There are 300 five-dimensional data points.
  x <- matrix(rnorm(kNumberOfDataPoints * kDimension, 0, 1), ncol = kDimension)

  # Randomly generate coefficients with different means.
  theta <- rbind(rnorm(kDimension, 0, 1), rnorm(kDimension, 2, 1))

  # Randomly generate response variables based on the segmented data and
  # corresponding coefficients
  y <- c(
    rbinom(kChangePointLocation, 1, 1 / (1 + exp(-x[1:kChangePointLocation, ] %*% theta[1, ]))),
    rbinom(kNumberOfDataPoints - kChangePointLocation, 1, 1 / (1 + exp(-x[(kChangePointLocation + 1):kNumberOfDataPoints, ] %*% theta[2, ])))
  )

  change_points_binomial_fastcpd <- suppressWarnings(fastcpd(
    formula = y ~ . - 1,
    data = data.frame(y = y, x = x),
    family = "binomial",
    segment_count = 5
  ))@cp_set

  change_points_binomial_fastcpd_sanity <- segd_binomial(
    cbind(x, y), (kDimension + 1) * log(kNumberOfDataPoints) / 2,
    B = 5
  )$cp

  expect_equal(change_points_binomial_fastcpd, change_points_binomial_fastcpd_sanity)

  change_points_binomial_fastcpd_vanilla_sanity <- pelt_vanilla_binomial(
    cbind(x, y), (kDimension + 1) * log(kNumberOfDataPoints) / 2
  )$cp
  expect_equal(change_points_binomial_fastcpd_vanilla_sanity, c(0, kChangePointLocation))
})
