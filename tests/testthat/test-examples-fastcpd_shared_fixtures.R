testthat::test_that(
  "shared deterministic fixtures have the same R-facing contract", {
    source("examples/fastcpd_shared_fixtures.R")

    for (case_id in names(shared_detector_results)) {
      testthat::expect_equal(
        shared_detector_results[[case_id]]@cp_set,
        shared_detector_expected_cp[[case_id]],
        info = case_id
      )
    }

    for (case_id in names(shared_variance_results)) {
      case <- shared_variance_cases[[case_id]]
      result <- shared_variance_results[[case_id]]
      expected <- as.numeric(case$expected_value)
      tolerance <- as.numeric(case$tolerance)

      if (case$operation == "estimate_variance_arma") {
        order <- parse_shared_order(case$order)
        testthat::expect_equal(
          nrow(result$direct$table),
          prod(order),
          info = case_id
        )
        testthat::expect_equal(
          rownames(result$direct$table),
          paste0("AR(", seq_len(prod(order)), ")"),
          info = case_id
        )
        direct_value <- result$direct$sigma2_bic
        generic_value <- result$generic$sigma2_bic
      } else {
        direct_value <- as.numeric(result$direct)
        generic_value <- as.numeric(result$generic)
      }

      testthat::expect_equal(
        direct_value,
        expected,
        tolerance = tolerance,
        info = case_id
      )
      testthat::expect_equal(
        generic_value,
        direct_value,
        tolerance = tolerance,
        info = case_id
      )
    }
  }
)
