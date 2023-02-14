#' Solve logistic/poisson regression using Gradient Descent Extension to the
#' multivariate case
#'
#' @param data A data frame containing the data to be segmented.
#' @param theta Estimate of the parameters. If null, the function will estimate
#'   the parameters.
#' @param family Family of the model.
#' @param lambda TODO
#'
#' @return Negative log likelihood of the corresponding data with the given
#'   family.
negative_log_likelihood <- function(data, theta, family, lambda = NULL) {
  data <- as.matrix(data)

  if (is.null(theta) && family == "gaussian") {
    # Estimate theta in gaussian family
    out <- glmnet::glmnet(
      as.matrix(data[, -1]), data[, 1],
      family = family, lambda = lambda
    )
    stats::deviance(out) / 2
  } else if (is.null(theta)) {
    # Estimate theta in binomial/poisson family
    out <- fastglm::fastglm(as.matrix(data[, -1]), data[, 1], family)
    out$deviance / 2
  } else if (family == "gaussian") {
    # Calculate negative log likelihood in gaussian family
    penalty <- lambda * sum(abs(theta))
    sum((data[, 1] - data[, -1, drop = FALSE] %*% theta)^2) / 2 + penalty
  } else if (family == "binomial") {
    # Calculate negative log likelihood in binomial family
    u <- c(data[, -1, drop = FALSE] %*% theta)
    sum(-data[, 1] * u + log(1 + exp(u)))
  } else {
    # Calculate negative log likelihood in poisson family
    u <- c(data[, -1, drop = FALSE] %*% theta)
    sum(-data[, 1] * u + exp(u) + lfactorial(data[, 1, drop = FALSE]))
  }
}
