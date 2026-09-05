testthat::test_that(
  "shared deterministic fixtures define the portable numerical contract", {
    source("examples/fastcpd_shared_fixtures.R")

    expect_shared_absolute_equal <- function(
      actual, expected, tolerance, info
    ) {
      actual <- as.numeric(actual)
      expected <- as.numeric(expected)
      testthat::expect_identical(
        is.na(actual),
        is.na(expected),
        info = paste(info, "missing-value positions")
      )
      comparable <- !is.na(actual) & !is.na(expected)
      if (!any(comparable)) return(invisible(actual))

      differences <- ifelse(
        actual[comparable] == expected[comparable],
        0,
        abs(actual[comparable] - expected[comparable])
      )
      worst <- which.max(differences)
      testthat::expect_true(
        differences[[worst]] <= tolerance,
        info = paste0(
          info,
          "\nactual: ", format(actual[comparable][[worst]], digits = 17),
          "\nexpected: ", format(expected[comparable][[worst]], digits = 17),
          "\nabsolute difference: ",
          format(differences[[worst]], digits = 17),
          "\nabsolute tolerance: ", format(tolerance, digits = 17)
        )
      )
      invisible(actual)
    }

    expected_families <- c(
      "mean", "variance", "meanvariance", "exponential", "lm", "lasso",
      "binomial", "poisson", "quantile", "var", "rank", "kcp", "ar",
      "arma", "arima", "garch"
    )
    testthat::expect_setequal(
      unique(vapply(shared_detector_cases, `[[`, character(1), "family")),
      expected_families
    )
    testthat::expect_setequal(
      unique(vapply(shared_confidence_cases, `[[`, character(1), "operation")),
      c("confint_bootstrap", "confint_profile", "confint_wald")
    )
    testthat::expect_s3_class(shared_var_p_response_error, "error")
    testthat::expect_match(
      conditionMessage(shared_var_p_response_error),
      "raw multivariate series"
    )
    testthat::expect_s3_class(shared_multivariate_wald_error, "error")
    testthat::expect_match(
      conditionMessage(shared_multivariate_wald_error),
      "multivariate LM"
    )
    testthat::expect_true(all(is.na(shared_separated_binomial_wald$se)))
    testthat::expect_true(all(is.na(shared_separated_binomial_wald$lower)))
    testthat::expect_true(all(is.na(shared_separated_binomial_wald$upper)))

    for (case_id in names(shared_detector_results)) {
      testthat::expect_equal(
        shared_detector_results[[case_id]]@cp_set,
        shared_detector_expected_cp[[case_id]],
        info = case_id
      )
    }

    for (index in seq_len(nrow(shared_expected_output_rows))) {
      row <- shared_expected_output_rows[index, , drop = FALSE]
      actual <- shared_numeric_outputs[[row$case_id]][[row$field]]
      expected <- parse_shared_expected_output(row)
      testthat::expect_identical(dim(actual), dim(expected), info = paste(
        row$case_id, row$field, "shape"
      ))
      expect_shared_absolute_equal(
        as.numeric(actual),
        as.numeric(expected),
        as.numeric(row$tolerance),
        paste(row$case_id, row$field)
      )
    }

    listed_files <- unique(vapply(
      shared_fixture_cases,
      function(case) case$data_file,
      character(1)
    ))
    listed_files <- listed_files[listed_files != "-"]
    testthat::expect_setequal(
      listed_files,
      basename(list.files(.shared_fixture_root, pattern = "[.]csv$"))
    )

    for (case_id in names(shared_variance_results)) {
      case <- shared_variance_cases[[case_id]]
      result <- shared_variance_results[[case_id]]
      expected <- as.numeric(case$expected_value)
      tolerance <- as.numeric(case$tolerance)

      if (case$operation == "estimate_variance_arma") {
        order <- parse_shared_numbers(case$order)
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

      expect_shared_absolute_equal(
        direct_value,
        expected,
        tolerance,
        case_id
      )
      expect_shared_absolute_equal(
        generic_value,
        direct_value,
        tolerance,
        case_id
      )
    }
  }
)
