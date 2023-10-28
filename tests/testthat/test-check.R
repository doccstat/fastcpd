testthat::test_that(
  "fastcpd.ts/fastcpd_ts", {
    testthat::expect_error(
      fastcpd.ts(),
      r"[The family should be one of "ar", "var", "ma", "arima", "garch".]"
    )

    testthat::expect_error(
      fastcpd.ts(family = "at"),
      r"[The family should be one of "ar", "var", "ma", "arima", "garch".
The provided family is "at".]"
    )

    testthat::expect_error(
      fastcpd_ts(family = "ar", order = c(0)),
      "The order should be a positive integer for AR family."
    )

    testthat::expect_error(
      fastcpd.ts(family = "ar", order = c(-1, 0, 0)),
      paste0(
        "The first element of the order should be a positive integer ",
        "for AR family."
      )
    )

    testthat::expect_error(
      fastcpd_ts(family = "arima", order = 1.1),
      "The order should be specified as a vector of length 3 for ARIMA family."
    )

    testthat::expect_error(
      fastcpd_ts(
        data = matrix(NA, 1, 2),
        family = "ar",
        order = 1
      ),
      "Data should be a univariate time series."
    )
  }
)

testthat::test_that(
  "invalid family provided", {
    testthat::expect_error(
      fastcpd(
        data = data.frame(y = 0, x = 0),
        family = "bin0mial"
      ),
      paste0(
        "The family should be one of \"lm\", \"binomial\", \"poisson\", ",
        "\"lasso\", \"ar\", \"var\", \"ma\", \"arima\", \"garch\", \"custom\".",
        "\nThe provided family is \"bin0mial\"."
      )
    )
  }
)

testthat::test_that(
  "custom family but no cost function provided", {
    testthat::expect_error(
      fastcpd(
        data = data.frame(y = 0, x = 0)
      ),
      "Please specify the cost function."
    )
  }
)

testthat::test_that(
  "linear regression without change points", {
    result <- fastcpd(
      formula = y ~ . - 1,
      data = data.frame(y = seq_len(100), x = seq_len(100)),
      family = "lm"
    )

    testthat::expect_length(result@cp_set, 0)
  }
)

testthat::test_that(
  "invalid order for ar family", {
    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "ar"
      ),
      paste0(
        "The first element of the order should be a positive integer ",
        "for AR family."
      )
    )

    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "ar",
        order = -0.1
      ),
      "The order should be a positive integer for AR family."
    )

    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "var",
        order = -0.1
      ),
      "The order should be a positive integer for VAR family."
    )

    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "arima",
        order = c(1, 0)
      ),
      "The order should be specified as a vector of length 3."
    )

    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "arima",
        order = c(0, 0, 0)
      ),
      "The order should have at least one non-zero element for ARIMA family."
    )

    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "arima",
        order = c(-0.1, 0, 0)
      ),
      "The order should be positive integers for ARIMA family."
    )
  }
)
