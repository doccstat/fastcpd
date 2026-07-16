testthat::test_that(
  "examples/fastcpd_arima.txt", {
    testthat::skip_if_not_installed("ggplot2")

    examples_arima <- readLines("examples/fastcpd_arima.txt")
    source(textConnection(paste(
      examples_arima[seq_len(length(examples_arima) - 2) + 1],
      collapse = "\n"
    )))

    testthat::expect_equal(result@cp_set, 41)
    mbic_adjustment <- log(41) / 2
    expected_costs <- c(
      20 * (log(2 * pi) + log(0.01) + 1) + mbic_adjustment,
      20 * (log(2 * pi) + log(4) + 1) + mbic_adjustment
    )
    testthat::expect_equal(result@cost_values, expected_costs)
    testthat::expect_equal(
      unname(as.matrix(result@thetas)),
      matrix(c(0.01, 4), 1)
    )
    testthat::expect_equal(which(is.na(result@residuals)), c(1, 42))
    testthat::expect_equal(
      c(result@residuals[!is.na(result@residuals)]),
      c(small_increments, large_increments[-1])
    )
    testthat::expect_error(
      detect_arima(x, c(0, 1, 0), include.mean = TRUE),
      "include.mean = TRUE"
    )
    testthat::expect_identical(
      formals(detect_arima)[["include.mean"]],
      FALSE
    )
    testthat::expect_identical(
      formals(fastcpd_arima)[["include.mean"]],
      FALSE
    )
    testthat::expect_identical(
      formals(fastcpd.arima)[["include.mean"]],
      FALSE
    )
  }
)
