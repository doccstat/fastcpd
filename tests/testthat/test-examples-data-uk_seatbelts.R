testthat::test_that(
  "examples/data-uk_seatbelts.R", {
    source("examples/data-uk_seatbelts.R")
    testthat::expect_equal(result_ar@cp_set, c(71, 158))
    testthat::expect_equal(result_lm@cp_set, c(50, 153))
  }
)
