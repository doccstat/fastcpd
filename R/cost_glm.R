#' Logistic regression.
#'
#' @param data TODO
#' @param family TODO
#'
#' @return TODO
#' @export
cost_glm <- function(data, family) {
  data <- as.matrix(data)
  p <- dim(data)[2] - 1
  out <- fastglm::fastglm(as.matrix(data[, seq_len(p)]), data[, p + 1], family)
  return(out$deviance / 2)
}
