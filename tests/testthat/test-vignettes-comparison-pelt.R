# Everything in this script is provided as is. The purpose of this script is to
# do a sanity check on the C++ implementation of `fastcpd`.

testthat::skip("Skip due to time limit on CRAN.")

# nolint start: script provided as is

testthat::test_that("logistic regression", {
  set.seed(1)
  p <- 5
  x <- matrix(rnorm(300 * p, 0, 1), ncol = p)

  # Randomly generate coefficients with different means.
  theta <- rbind(rnorm(p, 0, 1), rnorm(p, 2, 1))

  # Randomly generate response variables based on the segmented data and
  # corresponding coefficients
  y <- c(
    rbinom(125, 1, 1 / (1 + exp(-x[1:125, ] %*% theta[1, ]))),
    rbinom(300 - 125, 1, 1 / (1 + exp(-x[(125 + 1):300, ] %*% theta[2, ])))
  )

  change_points_binomial_fastcpd <- suppressWarnings(fastcpd.binomial(
    cbind(y, x),
    segment_count = 5,
    beta = "BIC",
    cost_adjustment = "BIC"
  ))@cp_set

  testthat::expect_equal(change_points_binomial_fastcpd, 125)

  warning_messages <- testthat::capture_warnings(
    change_points_binomial_fastcpd_vanilla <- fastcpd.binomial(
      cbind(y, x),
      segment_count = 5,
      vanilla_percentage = 1,
      beta = "BIC",
      cost_adjustment = "BIC"
    )@cp_set
  )

  testthat::expect_equal(
    sort(warning_messages),
    rep(c(
      "fit_glm: algorithm did not converge",
      "fit_glm: fitted probabilities numerically 0 or 1 occurred"
    ), c(6, 3984))
  )

  testthat::expect_equal(change_points_binomial_fastcpd_vanilla, 125)
})

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

  testthat::expect_equal(
    warning_messages, rep("fit_glm: fitted rates numerically 0 occurred", 1539)
  )

  testthat::expect_equal(
    change_points_poisson_fastcpd_vanilla,
    c(374, 752, 1133)
  )
})

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
    cost_adjustment = "BIC"
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
    cost_adjustment = "BIC"
  )@cp_set

  testthat::expect_equal(
    change_points_lasso_fastcpd_vanilla,
    c(103, 299, 510, 800, 900)
  )
})

# nolint end
