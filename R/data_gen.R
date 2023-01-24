#' Generate data from logistic regression models with change-points
#'
#' @param n TODO
#' @param d TODO
#' @param true.coef TODO
#' @param true.cp.loc TODO
#' @param Sigma TODO
#' @param family TODO
#' @param evar TODO
#'
#' @return TODO
#' @export
data_gen <- function(n, d, true.coef, true.cp.loc, Sigma, family, evar = NULL) {
  loc <- unique(c(0, true.cp.loc, n))
  if (dim(true.coef)[2] != length(loc) - 1) {
    stop("true.coef and true.cp.loc do not match")
  }
  x <- mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = Sigma)
  y <- NULL
  for (i in 1:(length(loc) - 1)) {
    if (family == "binomial") {
      Xb <- x[(loc[i] + 1):loc[i + 1], , drop = FALSE] %*% true.coef[, i, drop = FALSE]
      prob <- 1 / (1 + exp(-Xb))
      group <- stats::rbinom(length(prob), 1, prob)
      y <- c(y, group)
    } else if (family == "poisson") {
      mu <- exp(x[(loc[i] + 1):loc[i + 1], , drop = FALSE] %*% true.coef[, i, drop = FALSE])
      group <- stats::rpois(length(mu), mu)
      y <- c(y, group)
    } else if (family == "gaussian") {
      stopifnot(!is.null(evar))
      Xb <- x[(loc[i] + 1):loc[i + 1], , drop = FALSE] %*% true.coef[, i, drop = FALSE]
      add <- Xb + stats::rnorm(length(Xb), sd = sqrt(evar))
      y <- c(y, add)
    }
  }
  data <- cbind(x, y)
  true_cluster <- rep(1:(length(loc) - 1), diff(loc))
  result <- list(data, true_cluster)
  return(result)
}
