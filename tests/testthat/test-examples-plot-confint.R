testthat::test_that(
  "examples/plot-confint.R", {
    testthat::skip_if_not_installed("ggplot2")

    testthat::expect_no_error(source("examples/plot-confint.R"))

    testthat::expect_s3_class(cp_profile_interval, "fastcpd_confint")
    testthat::expect_s3_class(theta_wald_interval, "fastcpd_confint")

    # Plotting change-point intervals needs the fitted object carried on the
    # confint result; a plain data frame copy must fail loudly.
    stripped <- cp_profile_interval
    attr(stripped, "object") <- NULL
    testthat::expect_error(
      plot(stripped),
      "requires the fastcpd object"
    )
  }
)
