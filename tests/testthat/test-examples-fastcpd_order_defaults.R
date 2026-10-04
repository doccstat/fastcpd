testthat::test_that("omitted orders agree across public entry points", {
  source("examples/fastcpd_order_defaults.R")
  for (family in names(default_order_results)) {
    results <- default_order_results[[family]]
    for (entry in c("named", "generic", "explicit")) {
      testthat::expect_equal(
        results[[entry]]@order, default_order_cases[[family]], info = family
      )
      testthat::expect_false(results[[entry]]@cp_only, info = family)
      for (field in c("cp_set", "raw_cp_set", "cost_values", "residuals", "thetas")) {
        testthat::expect_equal(
          methods::slot(results[[entry]], field),
          methods::slot(results$explicit, field), info = paste(family, field)
        )
      }
    }
    if (family != "kcp") {
      testthat::expect_true(inherits(results$invalid_generic, "error"), info = family)
      testthat::expect_true(inherits(results$invalid_named, "error"), info = family)
    }
  }
  testthat::expect_equal(default_order_mean@order, c(0, 0, 0))
})
