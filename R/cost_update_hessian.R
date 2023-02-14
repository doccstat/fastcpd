#' Function to calculate the Hessian matrix at the current data.
#'
#' @param data A data frame containing the data to be segmented.
#' @param theta Estimated theta from the previous iteration.
#' @param hessian Hessian matrix from the previous iteration.
#' @param family Family of the model.
#'
#' @return Hessian at the current data.
cost_update_hessian <- function(
  data,
  theta,
  hessian,
  family,
  G = NULL
) {
  new_data_x <- data[-1]
  if (family == "binomial") {
    prob <- 1 / (1 + exp(-new_data_x %*% theta))
    hessian + (new_data_x %o% new_data_x) * c((1 - prob) * prob)
  } else if (family == "poisson") {
    prob <- exp(new_data_x %*% theta)
    hessian + (new_data_x %o% new_data_x) * min(c(prob), G)
  }
}
