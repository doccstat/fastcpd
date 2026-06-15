testthat::test_that(
  "examples/fastcpd_exponential_1.R", {
    source("examples/fastcpd_exponential_1.R")

    testthat::expect_equal(result@cp_set, c(300, 700))
    testthat::expect_equal(
      unname(unlist(result@thetas)),
      c(1.004055, 0.09939726, 0.9088465),
      tolerance = 1e-6
    )
  }
)
