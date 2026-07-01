#!/usr/bin/env Rscript
# Code generation script: generates 17 consolidated family translation units
# (one per family) from pelt.tmpl and segd.tmpl by token substitution.
#
# Each TU contains all variant functions for that family, reducing header
# parsing (fastcpd_template.h + Armadillo + Rcpp) from 246x to 17x with
# zero impact on generated machine code — every specialization is still a
# distinct template instantiation with its own if-constexpr path.
#
# Usage (from package root):
#   Rscript tools/codegen.R <out_dir> <src_dir>
#
# out_dir : directory in which to write the generated .cc files (default: src)
# src_dir : directory containing codegen/pelt.tmpl and codegen/segd.tmpl
#           (default: src)

args    <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[1] else "src"
src_dir <- if (length(args) >= 2) args[2] else "src"

# Remove all previously generated fastcpd_*.cc files (never touches
# fastcpd_impl.cc, which is a static source file, not generated).
old_files <- list.files(out_dir, pattern = "^fastcpd_.+\\.cc$", full.names = TRUE)
old_files <- old_files[basename(old_files) != "fastcpd_impl.cc"]
if (length(old_files) > 0) file.remove(old_files)

pelt_tmpl <- readLines(file.path(src_dir, "codegen", "pelt.tmpl"))
segd_tmpl <- readLines(file.path(src_dir, "codegen", "segd.tmpl"))

# Both templates share the same structure:
#   lines 1-2  : #include directives  (vary only by {POLICY_HEADER})
#   line  3    : blank separator
#   lines 4+   : function body        (vary by all other tokens)
# We emit the header block once per file, then append each variant's body.
TMPL_HEADER_LINES <- 3L  # includes the blank separator

substitute_tokens <- function(lines, subs) {
  for (token in names(subs)) {
    lines <- gsub(token, subs[[token]], lines, fixed = TRUE)
  }
  lines
}

# ---------------------------------------------------------------------------
# Family tables  (name, policy_header, policy_class)
# ---------------------------------------------------------------------------
FAMILIES_PELT <- list(
  c("mean",         "mean",         "MeanFamily"),
  c("mgaussian",    "mgaussian",    "MgaussianFamily"),
  c("variance",     "variance",     "VarianceFamily"),
  c("meanvariance", "meanvariance", "MeanvarianceFamily"),
  c("garch",        "garch",        "GarchFamily"),
  c("exponential",  "exponential",  "ExponentialFamily")
)

# p == 1 specialisations: scalar fast paths via kNDims=1.
FAMILIES_PELT_1D <- list(
  c("mean_1d",         "mean",         "MeanFamily"),
  c("variance_1d",     "variance",     "VarianceFamily"),
  c("meanvariance_1d", "meanvariance", "MeanvarianceFamily")
)

FAMILIES_SEGD <- list(
  c("arma",     "arma",     "ArmaFamily"),
  c("binomial", "binomial", "BinomialFamily"),
  c("custom",   "custom",   "CustomFamily"),
  c("gaussian", "gaussian", "GaussianFamily"),
  c("lasso",    "lasso",    "LassoFamily"),
  c("ma",       "ma",       "MaFamily"),
  c("poisson",  "poisson",  "PoissonFamily"),
  c("quantile", "quantile", "QuantileFamily")
)

# ---------------------------------------------------------------------------
# Variant axes  (cpp_value, suffix)
# ---------------------------------------------------------------------------
COST_ADJUSTMENTS     <- list(c("kBIC", ""), c("kMBIC", "_mbic"), c("kMDL", "_mdl"))
PROGRESS_VARIANTS    <- list(c("false", ""), c("true", "_prog"))
VANILLA_VARIANTS     <- list(c("false", ""), c("true", "_van"))
LINE_SEARCH_VARIANTS <- list(c("false", ""), c("true", "_ls"))

# ---------------------------------------------------------------------------
# Emit one consolidated .cc per PELT family.
# Contains 6 functions: 3 cost x 2 progress.
# ---------------------------------------------------------------------------
generate_pelt <- function(fam, k_n_dims) {
  tmpl_header <- pelt_tmpl[seq_len(TMPL_HEADER_LINES)]
  tmpl_body   <- pelt_tmpl[(TMPL_HEADER_LINES + 1L):length(pelt_tmpl)]

  out_lines <- substitute_tokens(tmpl_header, list("{POLICY_HEADER}" = fam[[2]]))

  for (cost in COST_ADJUSTMENTS) {
    for (prog in PROGRESS_VARIANTS) {
      func_suffix <- paste0(cost[[2]], prog[[2]])
      body <- substitute_tokens(tmpl_body, list(
        "{FAMILY_NAME}"       = fam[[1]],
        "{POLICY_CLASS}"      = fam[[3]],
        "{K_R_PROGRESS}"      = prog[[1]],
        "{K_VANILLA_ONLY}"    = "false",
        "{K_COST_ADJUSTMENT}" = cost[[1]],
        "{K_LINE_SEARCH}"     = "false",
        "{K_N_DIMS}"          = k_n_dims,
        "{FUNC_SUFFIX}"       = func_suffix
      ))
      out_lines <- c(out_lines, body, "")
    }
  }

  writeLines(out_lines, file.path(out_dir, paste0("fastcpd_", fam[[1]], ".cc")))
}

# ---------------------------------------------------------------------------
# Emit one consolidated .cc per SEGD family.
# Contains 24 functions: 3 cost x 2 vanilla x 2 line_search x 2 progress.
# ---------------------------------------------------------------------------
generate_segd <- function(fam) {
  tmpl_header <- segd_tmpl[seq_len(TMPL_HEADER_LINES)]
  tmpl_body   <- segd_tmpl[(TMPL_HEADER_LINES + 1L):length(segd_tmpl)]

  out_lines <- substitute_tokens(tmpl_header, list("{POLICY_HEADER}" = fam[[2]]))

  for (cost in COST_ADJUSTMENTS) {
    for (van in VANILLA_VARIANTS) {
      for (ls in LINE_SEARCH_VARIANTS) {
        for (prog in PROGRESS_VARIANTS) {
          func_suffix <- paste0(cost[[2]], van[[2]], ls[[2]], prog[[2]])
          body <- substitute_tokens(tmpl_body, list(
            "{FAMILY_NAME}"       = fam[[1]],
            "{POLICY_CLASS}"      = fam[[3]],
            "{K_R_PROGRESS}"      = prog[[1]],
            "{K_VANILLA_ONLY}"    = van[[1]],
            "{K_COST_ADJUSTMENT}" = cost[[1]],
            "{K_LINE_SEARCH}"     = ls[[1]],
            "{K_N_DIMS}"          = "-1",
            "{FUNC_SUFFIX}"       = func_suffix
          ))
          out_lines <- c(out_lines, body, "")
        }
      }
    }
  }

  writeLines(out_lines, file.path(out_dir, paste0("fastcpd_", fam[[1]], ".cc")))
}

# ---------------------------------------------------------------------------
# Generate: 9 PELT files + 8 SEGD files = 17 files total
# ---------------------------------------------------------------------------
for (fam in FAMILIES_PELT)    generate_pelt(fam, "-1")
for (fam in FAMILIES_PELT_1D) generate_pelt(fam, "1")
for (fam in FAMILIES_SEGD)    generate_segd(fam)

cat(sprintf("Generated 17 .cc files in %s\n", out_dir))
