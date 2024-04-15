require_namespace <- function(package_name) {
  requireNamespace(package_name, quietly = TRUE)
}

check_family <- function(family, allowed_families) {
  error_message <- paste0(
    "The family should be one of ",
    paste("\"", allowed_families, "\"", collapse = ", ", sep = ""),
    "."
  )
  if (is.null(family)) {
    stop(error_message)
  } else if (!(family %in% allowed_families)) {
    stop(paste0(error_message, "\nThe provided family is \"", family, "\"."))
  }
}

check_order <- function(order, family) {
  if (is.null(order)) {
    stop("Please refer to the documentation for the order of the model.")
  }
  switch(family,
    ar = check_ar_order(order),
    var = check_var_order(order),
    arima = check_arima_order(order),
    arma = check_arima_order(c(order[1], 0, order[2])),
    garch = check_garch_order(order)
  )
}

check_ar_order <- function(order) {  # nolint: cyclomatic complexity
  error_message <- c("The", "", "for AR family.")
  if (length(order) != 1 && length(order) != 3) {
    error_message[2] <- "order should be specified as a vector of length 1 or 3"
    stop(paste(error_message, collapse = " "))
  } else if (length(order) == 1 && (order <= 0 || order != floor(order))) {
    error_message[2] <- "order should be a positive integer"
    stop(paste(error_message, collapse = " "))
  } else if (
    length(order) == 3 && (order[1] <= 0 || order[1] != floor(order[1]))
  ) {
    error_message[2] <-
      "first element of the order should be a positive integer"
    stop(paste(error_message, collapse = " "))
  } else if (length(order) == 3 && (order[2] != 0 || order[3] != 0)) {
    error_message[2] <- "second and third elements of the order should be 0"
    stop(paste(error_message, collapse = " "))
  } else {
    TRUE
  }
}

check_var_order <- function(order) {
  error_message <- c("The", "", "for VAR family.")
  if (length(order) != 1) {
    error_message[2] <- "order should be specified as a single integer"
    stop(paste(error_message, collapse = " "))
  } else if (order <= 0 || order != floor(order)) {
    error_message[2] <- "order should be a positive integer"
    stop(paste(error_message, collapse = " "))
  } else {
    TRUE
  }
}

check_arima_order <- function(order) {
  error_message <- c("The", "", "for ARIMA family.")
  if (length(order) != 3) {
    error_message[2] <- "order should be specified as a vector of length 3"
    stop(paste(error_message, collapse = " "))
  } else if (any(order < 0) || any(order != floor(order))) {
    error_message[2] <- "order should be positive integers"
    stop(paste(error_message, collapse = " "))
  } else if (all(order == 0)) {
    error_message[2] <- "order should have at least one non-zero element"
    stop(paste(error_message, collapse = " "))
  } else {
    TRUE
  }
}

check_garch_order <- function(order) {
  error_message <- c("The", "", "for GARCH family.")
  if (length(order) != 2) {
    error_message[2] <- "order should be specified as a vector of length 2"
    stop(paste(error_message, collapse = " "))
  } else if (any(order < 0) || any(order != floor(order))) {
    error_message[2] <- "order should be positive integers"
    stop(paste(error_message, collapse = " "))
  } else if (all(order == 0)) {
    error_message[2] <- "order should have at least one non-zero element"
    stop(paste(error_message, collapse = " "))
  } else {
    TRUE
  }
}

check_cost <- function(cost, cost_gradient, cost_hessian, family) {  # nolint
  error_message <- c("Please", "")
  if (
    !(is.null(cost) && is.null(cost_gradient) && is.null(cost_hessian)) &&
      family != "custom"
  ) {
    error_message[2] <- "do not specify the cost function for built-in models."
    stop(paste(error_message, collapse = " "))
  } else if (family == "custom" && is.null(cost)) {
    error_message[2] <- "specify the cost function."
    stop(paste(error_message, collapse = " "))
  } else if (
    family == "custom" && is.null(cost_gradient) && !is.null(cost_hessian)
  ) {
    error_message[2] <-
      "specify the gradient function if the Hessian function is available."
    stop(paste(error_message, collapse = " "))
  } else if (
    family == "custom" && !is.null(cost_gradient) && is.null(cost_hessian)
  ) {
    error_message[2] <-
      "specify the Hessian function if the gradient function is available."
    stop(paste(error_message, collapse = " "))
  }
}

get_d <- function(data_, family) {
  if (family %in% c(
    "mean", "variance", "meanvariance", "ar", "arma", "arima", "garch", "var"
  )) {
    ncol(data_)
  } else {
    ncol(data_) - 1
  }
}

get_sigma_data <- function(data_, family, order, p, p_response) {  # nolint
  if (family == "ar") {
    y <- data_[p + seq_len(nrow(data_) - p), ]
    x <- matrix(NA, nrow(data_) - p, p)
    for (p_i in seq_len(p)) {
      x[, p_i] <- data_[(p - p_i) + seq_len(nrow(data_) - p), ]
    }
    data_ <- cbind(y, x)
  } else if (family == "var") {
    y <- data_[order + seq_len(nrow(data_) - order), ]
    x <- matrix(NA, nrow(data_) - order, order * ncol(data_))
    for (p_i in seq_len(order)) {
      x[, (p_i - 1) * ncol(data_) + seq_len(ncol(data_))] <-
        data_[(order - p_i) + seq_len(nrow(data_) - order), ]
    }
    data_ <- cbind(y, x)
  }

  sigma_ <- if (family == "mean") {
    variance.mean(data_)
  } else if (family == "var" || family == "lm" && p_response > 1) {
    as.matrix(Matrix::nearPD(variance.lm(data_, p_response))$mat)
  } else if (family == "lm" || family == "ar") {
    as.matrix(variance.lm(data_))
  } else {
    diag(1)
  }

  if (family == "mean") {
    sigma_inv <- solve(sigma_)
    chol_upper <- chol(sigma_inv)
    data1 <- tcrossprod(data_, chol_upper)
    data_ <- apply(cbind(data1, rowSums(data1^2)), 2, cumsum)
  } else if (family == "variance") {
    data_ <- data_ - colMeans(data_)
    data_ <- apply(data_, 1, tcrossprod)
    if (is.null(dim(data_)) || ncol(data_) == 1) {
      data_ <- matrix(data_)
    } else {
      data_ <- t(data_)
    }
    data_ <- apply(data_, 2, cumsum)
  } else if (family == "meanvariance") {
    data2 <- apply(data_, 1, tcrossprod)
    if (is.null(dim(data_)) || ncol(data_) == 1) {
      data2 <- matrix(data2)
    } else {
      data2 <- t(data2)
    }
    data_ <- apply(cbind(data_, data2), 2, cumsum)
  }

  list(sigma = sigma_, data = data_)
}

get_fastcpd_family <- function(family, p_response) {
  if (family %in% c(
    "binomial", "poisson", "lasso",
    "mean", "variance", "meanvariance", "arma"
  )) {
    family
  } else if (family == "lm" && p_response == 1 || family == "ar") {
    "gaussian"
  } else if (family == "var" || family == "lm" && p_response > 1) {
    "mgaussian"
  } else {
    "custom"
  }
}

get_vanilla_percentage <- function(vanilla_percentage, cost, fastcpd_family) {
  if (!is.null(cost) && length(formals(cost)) == 1 || fastcpd_family %in% c(
    "mean", "variance", "meanvariance", "arima", "garch", "mgaussian"
  )) {
    1
  } else {
    vanilla_percentage
  }
}

get_beta <- function(beta, p, n, fastcpd_family, sigma_) {
  if (is.character(beta)) {
    if (!(beta %in% c("BIC", "MBIC", "MDL"))) {
      stop("Invalid beta selection criterion provided.")
    }

    beta <- switch(
      beta,
      BIC = (p + 1) * log(n) / 2,
      MBIC = (p + 2) * log(n) / 2,
      MDL = (p + 2) * log2(n) / 2
    )

    # For linear regression models, an estimate of the variance is needed in the
    # cost function. The variance estimation is only for "lm" family with no
    # `beta` provided. Only estimate the variance for Gaussian family when
    # `beta` is null.
    if (fastcpd_family == "gaussian") {
      beta * c(sigma_)
    } else {
      beta
    }
  } else {
    beta
  }
}

get_pruning_coef <- function(
  pruning_coef_is_set,
  pruning_coef,
  cost_adjustment,
  fastcpd_family,
  p
) {
  if (!pruning_coef_is_set && (fastcpd_family %in% c("mgaussian", "lasso"))) {
    pruning_coef <- -Inf
  }
  if (!pruning_coef_is_set && cost_adjustment == "MBIC") {
    pruning_coef <- pruning_coef + p * log(2)
  } else if (!pruning_coef_is_set && cost_adjustment == "MDL") {
    pruning_coef <- pruning_coef + p * log2(2)
  }
  pruning_coef
}

get_p_response <- function(family, y, data) {
  if (family %in% c(
    "mean", "variance", "meanvariance", "arma", "arima", "garch"
  )) {
    0
  } else if (family == "var") {
    ncol(data)
  } else if (is.null(ncol(y))) {
    1
  } else {
    ncol(y)
  }
}

get_p <- function(data_, family, p_response, order, include_mean) {
  if (family == "mean") {
    ncol(data_)
  } else if (family == "variance") {
    ncol(data_)^2
  } else if (family == "meanvariance") {
    ncol(data_)^2 + ncol(data_)
  } else if (family == "ar") {
    order
  } else if (family == "arma") {
    sum(order) + 1
  } else if (family == "arima") {
    sum(order[-2]) + 1 + include_mean
  } else if (family == "garch") {
    sum(order) + 1
  } else if (family == "var") {
    order * p_response^2
  } else if (family == "lm" && p_response > 1) {
    (ncol(data_) - p_response) * p_response
  } else {
    ncol(data_) - 1
  }
}
