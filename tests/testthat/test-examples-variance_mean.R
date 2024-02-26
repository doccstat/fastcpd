testthat::test_that(
  "examples/variance_mean.R", {
    testthat::skip_if_not_installed("mvtnorm")

    source("examples/variance_mean.R")
    testthat::expect_equal(
      sigma,
      matrix(
        c(
          113.5442077,  -0.1433207,  -2.973472,
          -0.1433207, 106.1627730,   8.430343,
          -2.9734719,   8.4303433, 112.602989
        ), 3, 3, byrow = TRUE
      )
    )
  }
)
