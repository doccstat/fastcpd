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
  new_data_x <- data[-1]
  new_data_y <- data[1]
  if (family == "binomial") {
    c(-(new_data_y - 1 / (1 + exp(-new_data_x %*% theta)))) * new_data_x
  } else if (family == "poisson") {
    c(-(new_data_y - exp(new_data_x %*% theta))) * new_data_x
  }
}
