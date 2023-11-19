testthat::test_that(
  "examples/fastcpd_mean.R", {
    source("examples/fastcpd_mean.R")
    testthat::expect_equal(result@cp_set, c(300, 700))
  }
)
