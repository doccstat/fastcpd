testthat::test_that(
  "examples/interface_aliases.R", {
    source("examples/interface_aliases.R")

    testthat::expect_equal(
      detect_interface_result@cp_set,
      fastcpd_interface_result@cp_set
    )
    testthat::expect_equal(
      detect_generic_result@cp_set,
      fastcpd_interface_result@cp_set
    )
    testthat::expect_equal(
      estimate_variance_interface_result,
      variance_interface_result
    )
    testthat::expect_equal(
      estimate_variance_generic_result,
      variance_interface_result
    )

    detection_aliases <- c(
      detect = "fastcpd",
      detect_ar = "fastcpd_ar",
      detect_arima = "fastcpd_arima",
      detect_arma = "fastcpd_arma",
      detect_binomial = "fastcpd_binomial",
      detect_exponential = "fastcpd_exponential",
      detect_garch = "fastcpd_garch",
      detect_kcp = "fastcpd_kcp",
      detect_kernel = "fastcpd_kcp",
      detect_lasso = "fastcpd_lasso",
      detect_linear_regression = "fastcpd_lm",
      detect_lm = "fastcpd_lm",
      detect_logistic_regression = "fastcpd_binomial",
      detect_mean = "fastcpd_mean",
      detect_mean_variance = "fastcpd_meanvariance",
      detect_meanvariance = "fastcpd_meanvariance",
      detect_poisson = "fastcpd_poisson",
      detect_poisson_regression = "fastcpd_poisson",
      detect_quantile = "fastcpd_quantile",
      detect_quantile_regression = "fastcpd_quantile",
      detect_rank = "fastcpd_rank",
      detect_time_series = "fastcpd_ts",
      detect_ts = "fastcpd_ts",
      detect_var = "fastcpd_var",
      detect_variance = "fastcpd_variance"
    )
    for (alias in names(detection_aliases)) {
      testthat::expect_identical(
        getExportedValue("fastcpd", alias),
        getExportedValue("fastcpd", detection_aliases[[alias]])
      )
    }

    variance_aliases <- c(
      estimate_variance_arma = "variance_arma",
      estimate_variance_linear_regression = "variance_lm",
      estimate_variance_lm = "variance_lm",
      estimate_variance_mean = "variance_mean",
      estimate_variance_median = "variance_median"
    )
    for (alias in names(variance_aliases)) {
      testthat::expect_identical(
        getExportedValue("fastcpd", alias),
        getExportedValue("fastcpd", variance_aliases[[alias]])
      )
    }
  }
)
