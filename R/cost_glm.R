#' Logistic regression.
#'
#' @param data TODO
#' @param family TODO
#'
#' @return TODO
#' @export
cost_glm <- function(data, family = "binomial") {
  data <- as.matrix(data)
  p <- dim(data)[2] - 1
  x <- as.matrix(data[, 1:p])
  y <- data[, p + 1]
  out <- fastglm::fastglm(x = x, y = y, family = family)
  return(out$deviance / 2)
}
