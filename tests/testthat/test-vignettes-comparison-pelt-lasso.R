# Everything in this script is provided as is. The purpose of this script is to
# do a sanity check on the C++ implementation of `fastcpd`.

testthat::skip_on_cran()
testthat::skip_if_not_installed("mvtnorm")

# nolint start: script provided as is

# Generate data from penalized linear regression models with change-points
#' @param n Number of observations.
#' @param d Dimension of the covariates.
#' @param true.coef True regression coefficients.
#' @param true.cp.loc True change-point locations.
#' @param Sigma Covariance matrix of the covariates.
#' @param evar Error variance.
#' @keywords internal
#'
#' @noRd
#' @return A list containing the generated data and the true cluster
#'    assignments.
data_gen_lasso <- function(n, d, true.coef, true.cp.loc, Sigma, evar) {
  loc <- unique(c(0, true.cp.loc, n))
  if (dim(true.coef)[2] != length(loc) - 1) stop("true.coef and true.cp.loc do not match")
  x <- mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = Sigma)
  y <- NULL
  for (i in 1:(length(loc) - 1))
  {
    Xb <- x[(loc[i] + 1):loc[i + 1], , drop = FALSE] %*% true.coef[, i, drop = FALSE]
    add <- Xb + rnorm(length(Xb), sd = sqrt(evar))
    y <- c(y, add)
  }
  data <- cbind(x, y)
  true_cluster <- rep(1:(length(loc) - 1), diff(loc))
  result <- list(data, true_cluster)
  return(result)
}

testthat::test_that("penalized linear regression", {
  set.seed(1)
  n <- 1000
  s <- 3
  d <- 50
  evar <- 0.5
  Sigma <- diag(1, d)
  true.cp.loc <- c(100, 300, 500, 800, 900)
  seg <- length(true.cp.loc) + 1
  true.coef <- matrix(rnorm(seg * s), s, seg)
  true.coef <- rbind(true.coef, matrix(0, d - s, seg))
  out <- data_gen_lasso(n, d, true.coef, true.cp.loc, Sigma, evar)
  data <- out[[1]]
  beta <- log(n) / 2 # beta here has different meaning

  change_points_lasso_fastcpd <- fastcpd.lasso(
    cbind(data[, d + 1], data[, 1:d]),
    epsilon = 1e-5,
    beta = beta,
    cost_adjustment = "BIC",
    pruning_coef = 0
  )@cp_set

  testthat::expect_equal(
    change_points_lasso_fastcpd,
    c(100, 300, 520, 800, 901)
  )

  change_points_lasso_fastcpd_vanilla <- fastcpd.lasso(
    cbind(data[, d + 1], data[, 1:d]),
    vanilla_percentage = 1,
    epsilon = 1e-5,
    beta = beta,
    cost_adjustment = "BIC",
    pruning_coef = 0
  )@cp_set

  testthat::expect_equal(
    change_points_lasso_fastcpd_vanilla,
    c(103, 299, 510, 800, 900)
  )
})

# nolint end
