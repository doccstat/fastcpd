#' Cost function.
#'
#' @param data TODO
#' @param family TODO
#' @param lambda TODO
#'
#' @return TODO
#' @export
cost <- function(data, family, lambda = NULL) {
  data <- as.matrix(data)
  p <- dim(data)[2] - 1
  if (family == "gaussian") {
    stopifnot(!is.null(lambda))
    out <- glmnet::glmnet(
      as.matrix(data[, 1:p]), data[, p + 1], family = family, lambda = lambda
    )
    return(stats::deviance(out) / 2)
  } else if (family %in% c("binomial", "poisson")) {
    out <- fastglm::fastglm(as.matrix(data[, 1:p]), data[, p + 1], family)
    return(out$deviance / 2)
  } else {
    stop("family must be one of 'gaussian', 'binomial', or 'poisson'")
  }
}
