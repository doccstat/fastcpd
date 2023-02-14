#' Solve logistic/poisson regression using Gradient Descent Extension to the
#' multivariate case
#'
#' @param data TODO
#' @param b TODO
#' @param family TODO
#' @param lambda TODO
#'
#' @return TODO
#' @export
negative_log_likelihood <- function(data, b, family, lambda = NULL) {
  if (!(family %in% c("gaussian", "binomial", "poisson"))) {
    stop("family must be one of 'gaussian', 'binomial', or 'poisson'")
  }
  if (family == "gaussian" && is.null(lambda)) {
    stop("lambda must be specified for gaussian family")
  }
  data <- as.matrix(data)

  if (is.null(b) && family == "gaussian") {
    # Estimate b in gaussian family
    out <- glmnet::glmnet(
      as.matrix(data[, -1]), data[, 1],
      family = family, lambda = lambda
    )
    stats::deviance(out) / 2
  } else if (is.null(b)) {
    # Estimate b in binomial/poisson family
    out <- fastglm::fastglm(as.matrix(data[, -1]), data[, 1], family)
    out$deviance / 2
  } else if (family == "gaussian") {
    # Calculate negative log likelihood in gaussian family
    sum((data[, 1, drop = FALSE] - data[, -1, drop = FALSE] %*% b)**2 / 2 + lambda * abs(b))
  } else if (family == "binomial") {
    # Calculate negative log likelihood in binomial family
    u <- as.numeric(data[, -1, drop = FALSE] %*% b)
    sum(-data[, 1, drop = FALSE] * u + log(1 + exp(u)))
  } else {
    # Calculate negative log likelihood in poisson family
    u <- as.numeric(data[, -1, drop = FALSE] %*% b)
    sum(-data[, 1, drop = FALSE] * u + exp(u) + lfactorial(data[, 1, drop = FALSE]))
  }
}
