#' Function to calculate the Hessian matrix at the current data.
#'
#' @param data A data frame containing the data to be segmented.
#' @param theta Estimated theta from the previous iteration.
#' @param hessian Hessian matrix from the previous iteration.
#' @param family Family of the model.
#' @param min_prob Minimum probability to avoid numerical issues.
#'
#' @return Hessian at the current data.
cost_update_hessian <- function(
  data,
  theta,
  hessian,
  family,
  min_prob
) {
  data_x <- data[-1]
  if (family == "binomial") {

    prob <- 1 / (1 + exp(-data_x %*% theta))
    hessian + (data_x %o% data_x) * c((1 - prob) * prob)

  } else if (family == "poisson") {

    prob <- exp(data_x %*% theta)
    hessian + (data_x %o% data_x) * min(c(prob), min_prob)

  } else if (family %in% c("lasso", "gaussian")) {

    hessian + data_x %o% data_x

  } else {
    stop("family must be one of 'gaussian', 'binomial', or 'poisson'")
  }
}
