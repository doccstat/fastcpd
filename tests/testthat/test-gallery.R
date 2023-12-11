# Please make sure this file is in sync with "vignettes/gallery.Rmd"

testthat::test_that(
  "linear regression with multi-dimensional responses", {
    testthat::skip_if_not_installed("mvtnorm")
    set.seed(1)
    n <- 300
    p <- 3
    y_count <- 2
    x <- mvtnorm::rmvnorm(n, rep(0, p), diag(p))
    theta_0 <- array(NA, dim = c(3, y_count, 3))
    theta_0[, , 1] <- cbind(c(1, 1.2, -1), c(-1, 0, 0.5))
    theta_0[, , 2] <- cbind(c(-1, 0, 0.5), c(0.5, -0.3, 0.2))
    theta_0[, , 3] <- cbind(c(0.5, -0.3, 0.2), c(1, 1.2, -1))
    y <- rbind(
      x[1:100, ] %*% theta_0[, , 1],
      x[101:200, ] %*% theta_0[, , 2],
      x[201:n, ] %*% theta_0[, , 3]
    ) + matrix(rnorm(n * y_count), ncol = y_count)
    multi_response_linear_loss <- function(data) {
      x <- data[, (ncol(data) - p + 1):ncol(data)]
      y <- data[, 1:(ncol(data) - p)]

      if (nrow(data) <= p) {
        x_t_x <- diag(p)
      } else {
        x_t_x <- crossprod(x)
      }

      norm(y - x %*% solve(x_t_x, t(x)) %*% y, type = "F")^2 / 2
    }
    result <- fastcpd(
      formula = y ~ x - 1,
      data = data.frame(y = y, x = x),
      beta = (2 * p + 1) * log(n) / 2,
      cost = multi_response_linear_loss,
      cp_only = TRUE
    )

    testthat::expect_equal(result@cp_set, c(102, 195))
  }
)
