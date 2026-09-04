args <- commandArgs(trailingOnly = TRUE)
check_only <- identical(args, "--check")
if (length(args) > 1L || (length(args) == 1L && !check_only)) {
  stop("usage: Rscript generate_r_stochastic_fixtures.R [--check]")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
fixture_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1L]])))
} else {
  normalizePath("tests/fixtures")
}

write_fixture <- function(name, lines, directory = fixture_dir) {
  path <- file.path(directory, name)
  contents <- paste0(paste(lines, collapse = "\n"), "\n")
  if (check_only) {
    current <- paste0(paste(readLines(path, warn = FALSE), collapse = "\n"), "\n")
    if (!identical(current, contents)) {
      stop(name, " is not current; regenerate the R stochastic fixtures")
    }
  } else {
    writeChar(contents, path, eos = NULL)
  }
}

format_values <- function(values) {
  paste(sprintf("%.17g", values), collapse = ";")
}

rng_lines <- c(paste(
  "case_id", "seed", "operation", "population", "size", "replace",
  "expected", "tolerance", sep = "\t"
))
for (seed in c(1L, 7L, 123456789L)) {
  set.seed(seed)
  rng_lines <- c(rng_lines, paste(
    paste0("seed", seed, "_uniform"), seed, "uniform", "-", 6L, "-",
    format_values(stats::runif(6L)), 0, sep = "\t"
  ))
}
for (seed in c(1L, 7L, 123456789L)) {
  set.seed(seed)
  rng_lines <- c(rng_lines, paste(
    paste0("seed", seed, "_normal"), seed, "normal", "-", 6L, "-",
    format_values(stats::rnorm(6L)), "1e-14", sep = "\t"
  ))
}
for (seed in c(1L, 7L, 123456789L)) {
  set.seed(seed)
  rng_lines <- c(rng_lines, paste(
    paste0("seed", seed, "_replace"), seed, "sample", 8L, 12L, "true",
    paste(sample.int(8L, 12L, replace = TRUE), collapse = ";"), 0,
    sep = "\t"
  ))
}
for (seed in c(1L, 7L, 123456789L)) {
  set.seed(seed)
  rng_lines <- c(rng_lines, paste(
    paste0("seed", seed, "_no_replace"), seed, "sample", 20L, 10L, "false",
    paste(sample.int(20L, 10L), collapse = ";"), 0, sep = "\t"
  ))
}
write_fixture("r_rng.tsv", rng_lines)

kcp_input <- cbind(x = seq_len(24L) / 10, group = rep(c(-1, 1), 12L))
stochastic_dir <- file.path(fixture_dir, "stochastic")
if (!check_only) {
  dir.create(stochastic_dir, showWarnings = FALSE)
}
input_lines <- c(
  "x,group",
  apply(kcp_input, 1L, function(row) {
    paste(sprintf("%.17g", row), collapse = ",")
  })
)
write_fixture("kcp_seed_input.csv", input_lines, stochastic_dir)

set.seed(7L)
feature_count <- 8L
bandwidth <- 1.25
omega <- matrix(
  stats::rnorm(ncol(kcp_input) * feature_count, sd = 1 / bandwidth),
  nrow = ncol(kcp_input), ncol = feature_count
)
phase <- stats::runif(feature_count, 0, 2 * pi)
features <- sqrt(2 / feature_count) * cos(
  kcp_input %*% omega + matrix(
    phase, nrow(kcp_input), feature_count, byrow = TRUE
  )
)
feature_lines <- c(
  paste(paste0("feature_", seq_len(feature_count)), collapse = ","),
  apply(features, 1L, function(row) {
    paste(sprintf("%.17g", row), collapse = ",")
  })
)
write_fixture("kcp_seed_features.csv", feature_lines, stochastic_dir)

message(if (check_only) {
  "R stochastic fixtures are up to date."
} else {
  "Wrote R stochastic fixtures."
})
