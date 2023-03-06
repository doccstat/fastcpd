#' Function to calculate the gradient at the current data.
#'
#' @param data A data frame containing the data to be segmented.
#' @param theta Estimated theta from the previous iteration.
#' @param family Family of the model.
#'
#' @return Gradient at the current data.
cost_update_gradient <- function(
  data,
  theta,
  family
) {
  x <- data[-1]
  y <- data[1]
  if (family == "binomial") {
    c(-(y - 1 / (1 + exp(-x %*% theta)))) * x
  } else if (family == "poisson") {
    c(-(y - exp(x %*% theta))) * x
  } else if (family == "gaussian") {
    c(-(y - x %*% theta)) * x
  } else {
    stop("family must be one of 'gaussian', 'binomial', or 'poisson'")
  }
}
