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

#' Determine the number of arguments a `cost` function accepts.
#'
#' `cost` may be either an R closure -- where arity is `length(formals(cost))`
#' as before -- or a compiled cost function passed as an `externalptr` (built
#' with `Rcpp::XPtr`, see `?fastcpd`). External pointers carry no `formals`,
#' so compiled costs must instead carry an integer `fastcpd_cost_arity`
#' attribute: `1` for a PELT-style `cost(data)`, `2` for a SeGD-style
#' `cost(data, theta)`.
#' @noRd
cost_arity <- function(cost) {
  if (is.function(cost)) {
    length(formals(cost))
  } else if (inherits(cost, "externalptr")) {
    arity <- attr(cost, "fastcpd_cost_arity")
    if (is.null(arity) || !(arity %in% c(1L, 2L))) {
      stop(
        "Compiled (external pointer) cost functions must carry an integer ",
        "`fastcpd_cost_arity` attribute: 1 for a PELT-style cost(data), or ",
        "2 for a SeGD-style cost(data, theta). See `?fastcpd` for how to ",
        "build one with `Rcpp::XPtr`."
      )
    }
    as.integer(arity)
  } else {
    stop("`cost` must be a function or an external pointer (see `?fastcpd`).")
  }
}

get_pruning_coef <- function(
  pruning_coef_is_set,
  pruning_coef,
  cost_adjustment,
  fastcpd_family,
  p
) {
  # PELT pruning is only valid when segment costs satisfy a subadditivity
  # condition (Killick, Fearnhead & Eckley 2012, sec. 3.2): roughly, that
  # extending a segment cannot lower its optimal cost by more than a fixed
  # constant. This holds for the i.i.d.-type costs (mean/variance/regression)
  # PELT was designed for, but not for "mgaussian"/"lasso", nor for "garch",
  # whose likelihood is a recursive function of the whole segment: e.g. the
  # exact MLE NLL of a length-1 segment is 0 (no recursive terms to sum), so
  # short candidate segments can look spuriously cheap relative to the
  # threshold and cause the eventually-optimal split point to be discarded
  # before it is ever reconsidered. Disabling pruning keeps PELT exact
  # (falling back to optimal partitioning) for these families.
  if (!pruning_coef_is_set &&
      (fastcpd_family %in% c("mgaussian", "lasso", "garch"))) {
    pruning_coef <- -Inf
  }
  if (!pruning_coef_is_set && cost_adjustment == "MBIC") {
    pruning_coef <- pruning_coef + p * log(2)
  } else if (!pruning_coef_is_set && cost_adjustment == "MDL") {
    pruning_coef <- pruning_coef + p * log2(2)
  }
  pruning_coef
}

# Nearest positive definite matrix via eigendecomposition.
# Clips negative eigenvalues to .Machine$double.eps.
nearest_pd_ <- function(x) {
  eig <- eigen(x, symmetric = TRUE)
  eig$vectors %*%
    diag(pmax(eig$values, .Machine$double.eps), length(eig$values)) %*%
    t(eig$vectors)
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
