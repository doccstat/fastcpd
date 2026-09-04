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
  tolerance <- shared_fixture_cases[[case_id]]$tolerance
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
payload <- paste0(paste(lines, collapse = "\n"), "\n")

arguments <- commandArgs(TRUE)
check_only <- "--check" %in% arguments
if (check_only) {
  actual <- if (file.exists(output_path)) {
    paste0(paste(readLines(output_path, warn = FALSE), collapse = "\n"), "\n")
  } else {
    ""
  }
  if (!identical(actual, payload)) {
    cat("Out-of-date shared numerical outputs: expected_outputs.tsv\n")
    quit(status = 1L)
  }
  cat("Shared numerical outputs are up to date.\n")
} else {
  writeLines(lines, output_path, useBytes = TRUE)
  cat("Wrote shared numerical outputs to ", output_path, "\n", sep = "")
}
