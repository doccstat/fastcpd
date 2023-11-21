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
  } else {
    TRUE
  }
}

check_order <- function(order, family) {
  if (is.null(order)) {
    stop("Please refer to the documentation for the order of the model.")
  }
  switch(family,
    ar = check_ar_order(order),
    var = check_var_order(order),
    ma = check_ma_order(order),
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

check_ma_order <- function(order) {  # nolint: cyclomatic complexity
  error_message <- c("The", "", "for MA family.")
  if (length(order) != 1 && length(order) != 3) {
    error_message[2] <- "order should be specified as a vector of length 1 or 3"
    stop(paste(error_message, collapse = " "))
  } else if (length(order) == 1 && (order <= 0 || order != floor(order))) {
    error_message[2] <- "order should be a positive integer"
    stop(paste(error_message, collapse = " "))
  } else if (
    length(order) == 3 && (order[3] <= 0 || order[3] != floor(order[3]))
  ) {
    error_message[2] <-
      "third element of the order should be a positive integer"
    stop(paste(error_message, collapse = " "))
  } else if (length(order) == 3 && (order[1] != 0 || order[2] != 0)) {
    error_message[2] <- "first and second elements of the order should be 0"
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
  } else {
    TRUE
  }
}
