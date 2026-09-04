# Shared deterministic inputs and normalized outputs used by the R and Python
# suites. Resolve the path from both a source checkout and an installed package
# test directory.
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
  data <- as.matrix(read.csv(
    file.path(.shared_fixture_root, file),
    check.names = FALSE,
    stringsAsFactors = FALSE
  ))
  storage.mode(data) <- "double"
  data
}

.shared_manifest_columns <- c(
  "case_id", "data_file", "operation", "source_case", "family", "order",
  "beta", "cost_adjustment", "trim", "vanilla_percentage", "p_response",
  "variance_estimation", "random_state", "level", "B", "window",
  "expected_cp", "expected_value", "tolerance"
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

parse_shared_numbers <- function(value) {
  if (!nzchar(value) || identical(value, "-")) {
    return(numeric())
  }
  values <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
  if (anyNA(values)) stop("Invalid numeric fixture value: ", value)
  values
}

parse_shared_order <- function(value) {
  values <- parse_shared_numbers(value)
  if (length(values) == 1L) values[[1]] else values
}

parse_shared_beta <- function(value) {
  if (value %in% c("BIC", "MBIC", "MDL")) value else as.numeric(value)
}

parse_shared_change_points <- function(value) {
  if (!nzchar(value) || identical(value, "-")) numeric() else
    as.numeric(strsplit(value, ";", fixed = TRUE)[[1]])
}

shared_fixture_cases <- setNames(
  lapply(seq_len(nrow(shared_fixture_manifest)), function(index) {
    case <- as.list(shared_fixture_manifest[index, , drop = FALSE])
    case$data <- if (identical(case$data_file, "-")) NULL else
      read_shared_fixture(case$data_file)
    case
  }),
  shared_fixture_manifest$case_id
)

with_shared_seed <- function(seed, expression) {
  if (!length(seed)) return(force(expression))
  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed[[1]])
  force(expression)
}

run_shared_detector <- function(case) {
  arguments <- list(
    data = case$data,
    beta = parse_shared_beta(case$beta),
    cost_adjustment = case$cost_adjustment,
    trim = as.numeric(case$trim),
    vanilla_percentage = as.numeric(case$vanilla_percentage)
  )
  runner_name <- switch(
    case$family,
    mean = "detect_mean",
    variance = "detect_variance",
    meanvariance = "detect_meanvariance",
    exponential = "detect_exponential",
    lm = "detect_lm",
    lasso = "detect_lasso",
    binomial = "detect_binomial",
    poisson = "detect_poisson",
    quantile = "detect_quantile",
    var = "detect_var",
    rank = "detect_rank",
    kcp = "detect_kernel",
    ar = "detect_ar",
    arma = "detect_arma",
    arima = "detect_arima",
    garch = "detect_garch",
    stop("No R detector registered for family: ", case$family)
  )
  if (case$family %in% c("quantile", "var", "kcp", "ar", "arma", "arima", "garch")) {
    arguments$order <- parse_shared_order(case$order)
  }
  p_response <- parse_shared_numbers(case$p_response)
  if (length(p_response) && p_response[[1]] > 0) {
    arguments$p.response <- p_response[[1]]
  }
  variance_estimation <- parse_shared_numbers(case$variance_estimation)
  if (length(variance_estimation)) {
    arguments$variance_estimation <- if (length(variance_estimation) == 1L) {
      matrix(variance_estimation, 1L, 1L)
    } else {
      diag(variance_estimation)
    }
  }
  random_state <- parse_shared_numbers(case$random_state)
  result <- with_shared_seed(
    random_state,
    suppressWarnings(do.call(get(runner_name, mode = "function"), arguments))
  )
  # Store a stable public call so seeded bootstrap fixtures can replace the
  # data argument and refit through the same wrapper.
  result@call <- as.call(c(list(as.name(runner_name)), arguments))
  result
}

run_shared_variance_estimator <- function(case) {
  if (case$operation == "estimate_variance_arma") {
    order <- parse_shared_numbers(case$order)
    arguments <- list(data = case$data[, 1], p = order[1], q = order[2])
    return(list(
      direct = do.call(estimate_variance_arma, arguments),
      generic = do.call(
        estimate_variance,
        c(list(family = case$family), arguments)
      )
    ))
  }

  data <- if (case$family == "median") case$data[, 1] else case$data
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

normalize_shared_residuals <- function(result, case) {
  values <- as.numeric(result@residuals)
  if (!length(values)) return(matrix(numeric(), 0L, 0L))
  n <- nrow(result@data)
  order <- parse_shared_numbers(case$order)
  p_response <- parse_shared_numbers(case$p_response)
  response_count <- if (case$family == "var") {
    ncol(result@data)
  } else if (case$family == "lm" && length(p_response) && p_response[[1]] > 1) {
    p_response[[1]]
  } else if (case$family == "kcp") {
    as.integer(order[[1]])
  } else if (case$family %in% c("mean", "variance", "meanvariance")) {
    ncol(result@data)
  } else {
    1L
  }

  if (case$family == "var") {
    lag_count <- as.integer(order[[1]])
    native <- values[-seq_len(lag_count)]
    stopifnot(length(native) == (n - lag_count) * response_count)
    return(rbind(
      matrix(NA_real_, lag_count, response_count),
      matrix(native, n - lag_count, response_count)
    ))
  }
  stopifnot(length(values) == n * response_count)
  matrix(values, n, response_count)
}

run_shared_confidence <- function(case) {
  result <- shared_detector_results[[case$source_case]]
  arguments <- list(
    object = result,
    level = as.numeric(case$level)
  )
  if (case$operation == "confint_profile") {
    arguments$parm <- "cp"
    arguments$method <- "profile"
    arguments$window <- as.integer(case$window)
  } else if (case$operation == "confint_wald") {
    arguments$parm <- "theta"
    arguments$method <- "wald"
  } else if (case$operation == "confint_bootstrap") {
    arguments$parm <- "cp"
    arguments$method <- "bootstrap"
    arguments$B <- as.integer(case$B)
    arguments$seed <- as.integer(case$random_state)
  } else {
    stop("Unsupported confidence operation: ", case$operation)
  }
  suppressWarnings(do.call(stats::confint, arguments))
}

shared_confidence_cases <- Filter(
  function(case) startsWith(case$operation, "confint_"),
  shared_fixture_cases
)
shared_confidence_results <- lapply(
  shared_confidence_cases,
  run_shared_confidence
)

shared_numeric_outputs <- c(
  lapply(names(shared_detector_results), function(case_id) {
    result <- shared_detector_results[[case_id]]
    case <- shared_detector_cases[[case_id]]
    list(
      cp_set = as.numeric(result@cp_set),
      raw_cp_set = as.numeric(result@raw_cp_set),
      cost_values = as.numeric(result@cost_values),
      residuals = normalize_shared_residuals(result, case),
      thetas = as.matrix(result@thetas)
    )
  }) |> stats::setNames(names(shared_detector_results)),
  lapply(shared_confidence_results, function(result) {
    numeric_columns <- vapply(result, is.numeric, logical(1))
    lapply(result[numeric_columns], as.numeric)
  })
)

.shared_output_columns <- c(
  "case_id", "field", "shape", "values", "tolerance"
)
shared_expected_output_rows <- read.delim(
  file.path(.shared_fixture_root, "expected_outputs.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE,
  colClasses = "character",
  na.strings = NULL,
  quote = ""
)
if (!identical(names(shared_expected_output_rows), .shared_output_columns)) {
  stop("Unexpected shared numerical output columns")
}

parse_shared_expected_output <- function(row) {
  shape <- as.integer(strsplit(row$shape, ",", fixed = TRUE)[[1]])
  values <- if (identical(row$values, "-")) {
    numeric()
  } else {
    vapply(
      strsplit(row$values, ";", fixed = TRUE)[[1]],
      function(value) if (value %in% c("NA", "NaN")) NA_real_ else as.numeric(value),
      numeric(1)
    )
  }
  if (length(shape) == 2L) {
    values <- matrix(values, shape[1], shape[2], byrow = TRUE)
  }
  values
}
