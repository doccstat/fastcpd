testthat::test_that(
  "examples/interface_aliases.R", {
    source("examples/interface_aliases.R")

    testthat::expect_equal(
      detect_interface_result@cp_set,
      fastcpd_alias_result@cp_set
    )
    testthat::expect_equal(
      detect_generic_result@cp_set,
      detect_interface_result@cp_set
    )
    testthat::expect_equal(
      estimate_variance_interface_result,
      variance_alias_result
    )
    testthat::expect_equal(
      estimate_variance_generic_result,
      estimate_variance_interface_result
    )

    legacy_detection_aliases <- c(
      fastcpd = "detect",
      fastcpd_ar = "detect_ar",
      fastcpd.ar = "detect_ar",
      fastcpd_arima = "detect_arima",
      fastcpd.arima = "detect_arima",
      fastcpd_arma = "detect_arma",
      fastcpd.arma = "detect_arma",
      fastcpd_binomial = "detect_binomial",
      fastcpd.binomial = "detect_binomial",
      fastcpd_exponential = "detect_exponential",
      fastcpd.exponential = "detect_exponential",
      fastcpd_garch = "detect_garch",
      fastcpd.garch = "detect_garch",
      fastcpd_kcp = "detect_kernel",
      fastcpd.kcp = "detect_kernel",
      fastcpd_lasso = "detect_lasso",
      fastcpd.lasso = "detect_lasso",
      fastcpd_lm = "detect_lm",
      fastcpd.lm = "detect_lm",
      fastcpd_mean = "detect_mean",
      fastcpd.mean = "detect_mean",
      fastcpd_meanvariance = "detect_meanvariance",
      fastcpd.meanvariance = "detect_meanvariance",
      fastcpd_mv = "detect_meanvariance",
      fastcpd.mv = "detect_meanvariance",
      fastcpd_poisson = "detect_poisson",
      fastcpd.poisson = "detect_poisson",
      fastcpd_quantile = "detect_quantile",
      fastcpd.quantile = "detect_quantile",
      fastcpd_rank = "detect_rank",
      fastcpd.rank = "detect_rank",
      fastcpd_ts = "detect_time_series",
      fastcpd.ts = "detect_time_series",
      fastcpd_var = "detect_var",
      fastcpd.var = "detect_var",
      fastcpd_variance = "detect_variance",
      fastcpd.variance = "detect_variance"
    )
    for (alias in names(legacy_detection_aliases)) {
      testthat::expect_identical(
        getExportedValue("fastcpd", alias),
        getExportedValue("fastcpd", legacy_detection_aliases[[alias]])
      )
    }

    unified_detection_aliases <- c(
      detect_kcp = "detect_kernel",
      detect_linear_regression = "detect_lm",
      detect_logistic_regression = "detect_binomial",
      detect_mean_variance = "detect_meanvariance",
      detect_poisson_regression = "detect_poisson",
      detect_quantile_regression = "detect_quantile",
      detect_ts = "detect_time_series"
    )
    for (alias in names(unified_detection_aliases)) {
      testthat::expect_identical(
        getExportedValue("fastcpd", alias),
        getExportedValue("fastcpd", unified_detection_aliases[[alias]])
      )
    }

    legacy_variance_aliases <- c(
      variance_arma = "estimate_variance_arma",
      variance.arma = "estimate_variance_arma",
      variance_lm = "estimate_variance_linear_regression",
      variance.lm = "estimate_variance_linear_regression",
      variance_mean = "estimate_variance_mean",
      variance.mean = "estimate_variance_mean",
      variance_median = "estimate_variance_median",
      variance.median = "estimate_variance_median",
      estimate_variance_lm = "estimate_variance_linear_regression"
    )
    for (alias in names(legacy_variance_aliases)) {
      testthat::expect_identical(
        getExportedValue("fastcpd", alias),
        getExportedValue("fastcpd", legacy_variance_aliases[[alias]])
      )
    }
  }
)
