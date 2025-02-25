# Everything in this script is provided as is. The purpose of this script is to
# do a sanity check on the C++ implementation of `fastcpd`.

testthat::skip_on_cran()
testthat::skip_if_not_installed("mvtnorm")

# nolint start: script provided as is

# Generate data from poisson regression models with change-points
#' @param n Number of observations.
#' @param d Dimension of the covariates.
#' @param true.coef True regression coefficients.
#' @param true.cp.loc True change-point locations.
#' @param Sigma Covariance matrix of the covariates.
#' @keywords internal
#'
#' @noRd
#' @return A list containing the generated data and the true cluster
#'    assignments.
data_gen_poisson <- function(n, d, true.coef, true.cp.loc, Sigma) {
  loc <- unique(c(0, true.cp.loc, n))
  if (dim(true.coef)[2] != length(loc) - 1) stop("true.coef and true.cp.loc do not match")
  x <- mvtnorm::rmvnorm(n, mean = rep(0, d), sigma = Sigma)
  y <- NULL
  for (i in 1:(length(loc) - 1))
  {
    mu <- exp(x[(loc[i] + 1):loc[i + 1], , drop = FALSE] %*% true.coef[, i, drop = FALSE])
    group <- rpois(length(mu), mu)
    y <- c(y, group)
  }
  data <- cbind(x, y)
  true_cluster <- rep(1:(length(loc) - 1), diff(loc))
  result <- list(data, true_cluster)
  return(result)
}

testthat::test_that("poisson regression", {
  set.seed(1)
  n <- 1500
  d <- 5
  rho <- 0.9
  Sigma <- array(0, c(d, d))
  for (i in 1:d) {
    Sigma[i, ] <- rho^(abs(i - (1:d)))
  }
  delta <- c(5, 7, 9, 11, 13)
  a.sq <- 1
  delta.new <-
    delta * sqrt(a.sq) / sqrt(as.numeric(t(delta) %*% Sigma %*% delta))
  true.cp.loc <- c(375, 750, 1125)

  # regression coefficients
  true.coef <- matrix(0, nrow = d, ncol = length(true.cp.loc) + 1)
  true.coef[, 1] <- c(1, 1.2, -1, 0.5, -2)
  true.coef[, 2] <- true.coef[, 1] + delta.new
  true.coef[, 3] <- true.coef[, 1]
  true.coef[, 4] <- true.coef[, 3] - delta.new

  out <- data_gen_poisson(n, d, true.coef, true.cp.loc, Sigma)
  data <- out[[1]]
  g_tr <- out[[2]]
  beta <- log(n) * (d + 1) / 2

  change_points_poisson_fastcpd <- fastcpd.poisson(
    cbind(data[, d + 1], data[, 1:d]),
    beta = beta,
    cost_adjustment = "BIC",
    epsilon = 0.001,
    segment_count = 10
  )@cp_set

  testthat::expect_equal(
    change_points_poisson_fastcpd,
    c(380, 751, 1136, 1251)
  )

  warning_messages <- testthat::capture_warnings(
    change_points_poisson_fastcpd_vanilla <- fastcpd.poisson(
      cbind(data[, d + 1], data[, 1:d]),
      segment_count = 10,
      vanilla_percentage = 1,
      beta = beta,
      cost_adjustment = "BIC"
    )@cp_set
  )

  testthat::expect_length(warning_messages, 1766)

  testthat::expect_equal(
    change_points_poisson_fastcpd_vanilla,
    c(374, 752, 1133)
  )
})

# nolint end
