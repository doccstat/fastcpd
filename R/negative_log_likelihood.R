#' Solve logistic/poisson regression using Gradient Descent Extension to the
#' multivariate case
#'
#' @param data A data frame containing the data to be segmented.
#' @param theta Estimate of the parameters. If null, the function will estimate
#'   the parameters.
#' @param family Family of the model.
#' @param lambda Lambda for L1 regularization. Only used for lasso.
#'
#' @return Negative log likelihood of the corresponding data with the given
#'   family.
negative_log_likelihood <- function(data, theta, family, lambda, cv = FALSE) {
  data <- as.matrix(data)

  if (is.null(theta) && family == "lasso") {

    # Estimate theta in lasso family
    if (cv) {
      out <- glmnet::cv.glmnet(
        data[, -1, drop = FALSE],
        data[, 1],
        family = "gaussian"
      )
      return(list(theta = stats::coef(out, s = "lambda.1se")[-1], val = NULL))
    } else {
      out <- glmnet::glmnet(
        as.matrix(data[, -1]), data[, 1],
        family = "gaussian", lambda = lambda
      )
      return(list(theta = NULL, val = out$deviance / 2))
    }

  } else if (is.null(theta)) {

    # Estimate theta in binomial/poisson/gaussian family
    out <- fastglm::fastglm(as.matrix(data[, -1]), data[, 1], family)
    return(list(theta = out$coefficients, val = out$deviance / 2))

  } else if (family %in% c("lasso", "gaussian")) {

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
