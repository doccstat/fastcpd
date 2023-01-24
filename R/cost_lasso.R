#' Logistic regression.
#'
#' @param data TODO
#' @param lambda TODO
#' @param family TODO
#'
#' @return TODO
#' @export
cost_lasso <- function(data, lambda, family = "gaussian") {
  data <- as.matrix(data)
  p <- dim(data)[2] - 1
  out <- glmnet::glmnet(
    as.matrix(data[, 1:p]), data[, p + 1], family = family, lambda = lambda
  )
  return(stats::deviance(out) / 2)
}
