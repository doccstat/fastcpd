testthat::test_that(
  "the shared mean contract matches Python and C++", {
    source("examples/cross_language_mean.R")
    testthat::expect_equal(cross_language_mean_result@cp_set, 50)
  }
)
