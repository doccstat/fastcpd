testthat::test_that(
  "invalid family provided for time series data", {
    testthat::expect_error(
      fastcpd.ts(
        data = seq_len(10),
        family = "at",
        order = 1
      ),
      r"[
The family should be one of "ar", "var" or "arima",
while the provided family is at.]"
    )
  }
)

testthat::test_that(
  "invalid data provided for ar(1) model", {
    testthat::expect_error(
      fastcpd_ts(
        data = matrix(NA, 1, 2),
        family = "ar",
        order = 1
      ),
      "The data should be a univariate time series."
    )
  }
)

testthat::test_that(
  "ar(1) model using `fastcpd_ts`", {
    set.seed(1)
    n <- 1000
    p <- 1
    x <- rep(0, n + 1)
    for (i in 1:600) {
      x[i + 1] <- 0.6 * x[i] + rnorm(1)
    }
    for (i in 601:1000) {
      x[i + 1] <- 0.3 * x[i] + rnorm(1)
    }

    result <- fastcpd_ts(x, "ar", 1)

    testthat::expect_equal(result@cp_set, 609)
  }
)

testthat::test_that(
  "ma(4) model with `forecast::Arima` calculation of ML", {
    testthat::skip_if_not_installed("forecast")
    set.seed(1)
    n <- 400
    p <- 4
    time_series <- rep(0, n)

    lambda1 <- c(1, 3 / 2, 2i, -2i)
    param1 <- 1
    for (lam in lambda1) {
      param1 <- convolve(param1, c(1, -lam), type = "open")
    }
    ma_coef1 <- Re(param1[-1] / param1[1])

    testthat::expect_equal(ma_coef1, c(-5 / 3, 11 / 12, -5 / 12, 1 / 6))

    time_series[1:200] <-
      arima.sim(n = 200, model = list(ma = ma_coef1, sd = 0.1))

    lambda2 <- c(3, 1, -3i / 2, 3i / 2)
    param2 <- 1
    for (lam in lambda2) {
      param2 <- convolve(param2, c(1, -lam), type = "open")
    }
    ma_coef2 <- Re(param2[-1] / param2[1])

    testthat::expect_equal(ma_coef2, c(-4 / 3, 7 / 9, -16 / 27, 4 / 27))

    time_series[201:400] <-
      arima.sim(n = 200, model = list(ma = ma_coef2, sd = 0.3))

    result <- fastcpd.ts(
      data = data.frame(x = time_series),
      family = "arima",
      order = c(0, 0, 4)
    )

    testthat::expect_equal(result@cp_set, 200)
  }
)
