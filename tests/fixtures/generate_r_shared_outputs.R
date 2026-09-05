#!/usr/bin/env Rscript

# Generate the normalized numerical output contract from the R reference
# implementation. Python reads the committed TSV independently.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_argument)) stop("Unable to locate this script")
script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]))
fixture_root <- dirname(script_path)
package_root <- normalizePath(file.path(fixture_root, "..", ".."))
output_path <- file.path(fixture_root, "expected_outputs.tsv")

devtools::load_all(package_root, quiet = TRUE)
old_working_directory <- setwd(package_root)
on.exit(setwd(old_working_directory), add = TRUE)
source(file.path("tests", "testthat", "examples", "fastcpd_shared_fixtures.R"))

format_shared_value <- function(value) {
  if (is.na(value)) return("NA")
  if (is.nan(value)) return("NaN")
  if (is.infinite(value)) return(if (value > 0) "Inf" else "-Inf")
  sprintf("%.17g", value)
}

shared_output_tolerance <- function(case_id, field) {
  case <- shared_fixture_cases[[case_id]]
  if (
    identical(case_id, "binomial_wald") &&
      field %in% c("lower", "upper")
  ) {
    # Preserve the former combined relative/absolute Python allowance as one
    # explicit absolute bound for the accumulated IRLS/Wald rounding error.
    return("2.2e-7")
  }
  if (
    identical(case$operation, "confint_bootstrap") &&
      identical(field, "detection_rate") &&
      identical(shared_fixture_cases[[case$source_case]]$family, "kcp")
  ) {
    # A libm-sensitive KCP boundary can change one bootstrap detection across
    # operating systems. Include a few ulps so both adjacent B-count ratios
    # are inside the intended one-draw absolute tolerance.
    return(format_shared_value(
      (1 + 8 * .Machine$double.eps) / as.integer(case$B)
    ))
  }
  case$tolerance
}

serialize_shared_output <- function(case_id, field, value) {
  dimensions <- dim(value)
  if (is.null(dimensions)) {
    shape <- as.character(length(value))
    flattened <- as.numeric(value)
  } else {
    if (length(dimensions) != 2L) stop("Only vectors and matrices are supported")
    shape <- paste(dimensions, collapse = ",")
    flattened <- as.numeric(t(value))
  }
  values <- if (length(flattened)) {
    paste(vapply(flattened, format_shared_value, character(1)), collapse = ";")
  } else {
    "-"
  }
  tolerance <- shared_output_tolerance(case_id, field)
  paste(case_id, field, shape, values, tolerance, sep = "\t")
}

lines <- "case_id\tfield\tshape\tvalues\ttolerance"
for (case_id in names(shared_numeric_outputs)) {
  for (field in names(shared_numeric_outputs[[case_id]])) {
    lines <- c(lines, serialize_shared_output(
      case_id,
      field,
      shared_numeric_outputs[[case_id]][[field]]
    ))
  }
}
arguments <- commandArgs(TRUE)
check_only <- "--check" %in% arguments
if (check_only) {
  expected_keys <- paste(
    shared_expected_output_rows$case_id,
    shared_expected_output_rows$field,
    sep = "\t"
  )
  actual_keys <- unlist(lapply(names(shared_numeric_outputs), function(case_id) {
    paste(case_id, names(shared_numeric_outputs[[case_id]]), sep = "\t")
  }))
  stale <- if (!identical(expected_keys, actual_keys)) {
    "case/field inventory"
  } else {
    character()
  }

  if (!length(stale)) {
    for (index in seq_len(nrow(shared_expected_output_rows))) {
      row <- shared_expected_output_rows[index, , drop = FALSE]
      actual <- shared_numeric_outputs[[row$case_id]][[row$field]]
      expected <- parse_shared_expected_output(row)
      intended_tolerance <- as.numeric(shared_output_tolerance(
        row$case_id, row$field
      ))
      tolerance <- as.numeric(row$tolerance)
      label <- paste(row$case_id, row$field)

      if (!identical(dim(actual), dim(expected))) {
        stale <- c(stale, paste(label, "shape"))
        next
      }
      actual <- as.numeric(actual)
      expected <- as.numeric(expected)
      if (
        !identical(tolerance, intended_tolerance) ||
          !identical(is.na(actual), is.na(expected))
      ) {
        stale <- c(stale, label)
        next
      }
      comparable <- !is.na(actual) & !is.na(expected)
      differences <- ifelse(
        actual[comparable] == expected[comparable],
        0,
        abs(actual[comparable] - expected[comparable])
      )
      if (length(differences) && max(differences) > tolerance) {
        stale <- c(stale, label)
      }
    }
  }

  if (length(stale)) {
    cat(
      "Out-of-date shared numerical outputs: ",
      paste(stale, collapse = ", "),
      "\n",
      sep = ""
    )
    quit(status = 1L)
  }
  cat("Shared numerical outputs are up to date.\n")
} else {
  writeLines(lines, output_path, useBytes = TRUE)
  cat("Wrote shared numerical outputs to ", output_path, "\n", sep = "")
}
