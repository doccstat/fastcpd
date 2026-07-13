testthat::test_that(
  "examples/confint.R", {
    source("examples/confint.R")

    testthat::expect_gt(length(result@cp_set), 0)

    testthat::expect_s3_class(cp_profile_interval, "fastcpd_confint")
    testthat::expect_equal(cp_profile_interval$estimate, result@cp_set)
    testthat::expect_true(all(cp_profile_interval$lower <= result@cp_set))
    testthat::expect_true(all(cp_profile_interval$upper >= result@cp_set))

    testthat::expect_s3_class(theta_wald_interval, "fastcpd_confint")
    testthat::expect_equal(
      nrow(theta_wald_interval),
      ncol(result@thetas)
    )
    testthat::expect_true(all(
      theta_wald_interval$lower <= theta_wald_interval$estimate
    ))
    testthat::expect_true(all(
      theta_wald_interval$upper >= theta_wald_interval$estimate
    ))

    testthat::expect_s3_class(cp_bootstrap_interval, "fastcpd_confint")
    testthat::expect_equal(cp_bootstrap_interval$estimate, result@cp_set)
    testthat::expect_true(all(cp_bootstrap_interval$detection_rate >= 0))
    testthat::expect_true(all(cp_bootstrap_interval$detection_rate <= 1))
    testthat::expect_true(all(
      is.na(cp_bootstrap_interval$lower) |
        cp_bootstrap_interval$lower <= cp_bootstrap_interval$upper
    ))

    testthat::expect_error(
      confint(variance_result, parm = "theta", method = "wald"),
      "Wald intervals are not implemented"
    )
  }
)
