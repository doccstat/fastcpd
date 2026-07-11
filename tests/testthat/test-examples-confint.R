testthat::test_that(
  "examples/confint.R", {
    source("examples/confint.R")

    testthat::expect_length(result@cp_set, 1)

    testthat::expect_equal(cp_profile_interval$estimate, result@cp_set)
    testthat::expect_true(cp_profile_interval$lower <= result@cp_set)
    testthat::expect_true(cp_profile_interval$upper >= result@cp_set)

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
