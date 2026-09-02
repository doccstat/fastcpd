# Shared deterministic inputs used by the R and Python test suites.
# Resolve the path both from a source checkout (tests/testthat) and from an
# installed R package test directory.
.shared_fixture_candidates <- c(
  file.path("..", "fixtures"),
  file.path("tests", "fixtures")
)
.shared_fixture_root <- .shared_fixture_candidates[
  vapply(.shared_fixture_candidates, dir.exists, logical(1))
][1]
if (is.na(.shared_fixture_root)) {
  stop("Unable to locate tests/fixtures")
}

read_shared_fixture <- function(file) {
  read.csv(
    file.path(.shared_fixture_root, file),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

.shared_manifest_columns <- c(
  "case_id", "data_file", "operation", "family", "order", "beta",
  "cost_adjustment", "trim", "vanilla_percentage", "expected_cp",
  "expected_value", "tolerance"
)
shared_fixture_manifest <- read.delim(
  file.path(.shared_fixture_root, "manifest.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character",
  na.strings = NULL,
  quote = ""
)
if (!identical(names(shared_fixture_manifest), .shared_manifest_columns)) {
  stop("Unexpected shared fixture manifest columns")
}
if (anyDuplicated(shared_fixture_manifest$case_id)) {
  stop("Shared fixture case_id values must be unique")
}

parse_shared_order <- function(value) {
  if (!nzchar(value)) {
    return(integer())
  }
  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

parse_shared_beta <- function(value) {
  if (value %in% c("BIC", "MBIC", "MDL")) {
    return(value)
  }
  as.numeric(value)
}

parse_shared_change_points <- function(value) {
  if (!nzchar(value)) {
    return(numeric())
  }
  as.numeric(strsplit(value, ";", fixed = TRUE)[[1]])
}

shared_fixture_cases <- setNames(
  lapply(seq_len(nrow(shared_fixture_manifest)), function(index) {
    case <- as.list(shared_fixture_manifest[index, , drop = FALSE])
    case$data <- read_shared_fixture(case$data_file)
    case
  }),
  shared_fixture_manifest$case_id
)

run_shared_detector <- function(case) {
  arguments <- list(
    data = case$data,
    beta = parse_shared_beta(case$beta),
    cost_adjustment = case$cost_adjustment,
    trim = as.numeric(case$trim),
    vanilla_percentage = as.numeric(case$vanilla_percentage)
  )
  runner <- switch(
    case$family,
    mean = detect_mean,
    variance = detect_variance,
    meanvariance = detect_meanvariance,
    exponential = detect_exponential,
    lm = detect_lm,
    rank = detect_rank,
    arima = detect_arima,
    stop("No R detector registered for family: ", case$family)
  )
  if (case$family == "arima") {
    arguments$data <- case$data[[1]]
    arguments$order <- parse_shared_order(case$order)
  }
  if (case$family == "lm") {
    return(suppressWarnings(do.call(runner, arguments)))
  }
  do.call(runner, arguments)
}

run_shared_variance_estimator <- function(case) {
  if (case$operation == "estimate_variance_arma") {
    order <- parse_shared_order(case$order)
    arguments <- list(data = case$data[[1]], p = order[1], q = order[2])
    return(list(
      direct = do.call(estimate_variance_arma, arguments),
      generic = do.call(
        estimate_variance,
        c(list(family = case$family), arguments)
      )
    ))
  }

  data <- if (case$family == "median") case$data[[1]] else case$data
  runner <- switch(
    case$family,
    mean = estimate_variance_mean,
    median = estimate_variance_median,
    lm = estimate_variance_linear_regression,
    stop("No R variance estimator registered for family: ", case$family)
  )
  list(
    direct = runner(data),
    generic = estimate_variance(data, family = case$family)
  )
}

shared_detector_cases <- Filter(
  function(case) case$operation == "detect",
  shared_fixture_cases
)
shared_detector_results <- lapply(shared_detector_cases, run_shared_detector)
shared_detector_expected_cp <- lapply(
  shared_detector_cases,
  function(case) parse_shared_change_points(case$expected_cp)
)

shared_variance_cases <- Filter(
  function(case) startsWith(case$operation, "estimate_variance"),
  shared_fixture_cases
)
shared_variance_results <- lapply(
  shared_variance_cases,
  run_shared_variance_estimator
)
