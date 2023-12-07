testthat::test_that(
  "fastcpd.ts/fastcpd_ts", {
    testthat::expect_error(
      fastcpd.ts(),
      paste0(
        "The family should be one of \"ar\", \"var\", \"ma\", \"arima\", ",
        "\"arma\", \"garch\"."
      )
    )

    testthat::expect_error(
      fastcpd.ts(family = "at"),
      paste0(
        "The family should be one of \"ar\", \"var\", \"ma\", \"arima\", ",
        "\"arma\", \"garch\".\n",
        "The provided family is \"at\"."
      )
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
  "order should be specified as a vector of length 1 or 3", {
    testthat::expect_error(
      fastcpd.ts(seq_len(5), "ar", rep(2, 2)),
      paste0(
        "The order should be specified as a vector of length 1 or 3 ",
        "for AR family."
      )
    )
  }
)

testthat::test_that(
  "second and third elements of the order should be 0", {
    testthat::expect_error(
      fastcpd.ts(seq_len(5), "ar", c(1, 1, 0)),
      paste0(
        "The second and third elements of the order should be 0 ",
        "for AR family."
      )
    )
  }
)

testthat::test_that(
  "order should be specified as a single integer", {
    testthat::expect_error(
      fastcpd.ts(seq_len(5), "var", c(1, 1)),
      paste0(
        "The order should be specified as a single integer ",
        "for VAR family."
      )
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
        "\"lasso\", \"mlasso\", ",
        "\"mean\", \"variance\", \"meanvariance\", \"mv\", ",
        "\"arma\", \"ar\", \"var\", \"ma\", \"arima\", \"garch\", \"custom\".",
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
  "check_cost", {
    testthat::expect_error(
      fastcpd(
        data = data.frame(y = 0, x = 0),
        family = "lm",
        cost = function(data) NULL
      ),
      "Please do not specify the cost function for built-in models."
    )

    testthat::expect_error(
      fastcpd(
        data = data.frame(y = 0, x = 0),
        cost = function(data, theta) NULL,
        cost_gradient = function(data, theta) NULL
      ),
      paste0(
        "Please specify the Hessian function ",
        "if the gradient function is available."
      )
    )

    testthat::expect_error(
      fastcpd(
        data = data.frame(y = 0, x = 0),
        cost = function(data, theta) NULL,
        cost_hessian = function(data, theta) NULL
      ),
      paste0(
        "Please specify the gradient function ",
        "if the Hessian function is available."
      )
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
      fastcpd.ts(
        data = seq_len(5),
        family = "ar",
        order = NULL
      ),
      "Please refer to the documentation for the order of the model."
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
        family = "ma",
        order = c(1, 0)
      ),
      paste0(
        "The order should be specified as a vector of length 1 or 3 ",
        "for MA family."
      )
    )

    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "ma",
        order = 0
      ),
      "The order should be a positive integer for MA family."
    )

    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "ma",
        order = c(0, 0, -1)
      ),
      paste0(
        "The third element of the order should be a positive integer ",
        "for MA family."
      )
    )

    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "ma",
        order = c(1, 0, 1)
      ),
      "The first and second elements of the order should be 0 for MA family."
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

    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "garch",
        order = 1
      ),
      "The order should be specified as a vector of length 2 for GARCH family."
    )

    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "garch",
        order = c(0.1, 0.1)
      ),
      "The order should be positive integers for GARCH family."
    )

    testthat::expect_error(
      fastcpd(
        formula = ~ x - 1,
        data = data.frame(x = 0),
        family = "garch",
        order = c(0, 0)
      ),
      "The order should have at least one non-zero element for GARCH family."
    )
  }
)
