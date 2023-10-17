testthat::test_that(
  "invalid family provided for time series data", {
    testthat::expect_error(
      fastcpd.ts(
        data = seq_len(10),
        family = "at",
        order = 1
      ),
      r"[
The family should be one of "ar" or "var"
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
