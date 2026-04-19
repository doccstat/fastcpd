#!/usr/bin/env Rscript
# Code generation script: generates all 216 family translation units from
# pelt.tmpl and segd.tmpl by token substitution.
#
# Usage (from package root):
#   Rscript tools/codegen.R <out_dir> <src_dir>
#
# out_dir : directory in which to write the generated .cc files (default: src)
# src_dir : directory containing codegen/pelt.tmpl and codegen/segd.tmpl
#           (default: src)

args <- commandArgs(trailingOnly = TRUE)
out_dir  <- if (length(args) >= 1) args[1] else "src"
src_dir  <- if (length(args) >= 2) args[2] else "src"

pelt_tmpl <- readLines(file.path(src_dir, "codegen", "pelt.tmpl"))
segd_tmpl <- readLines(file.path(src_dir, "codegen", "segd.tmpl"))

generate <- function(tmpl, subs, outfile) {
  content <- tmpl
  for (token in names(subs)) {
    content <- gsub(token, subs[[token]], content, fixed = TRUE)
  }
  writeLines(content, outfile)
}

# ---------------------------------------------------------------------------
# Family tables  (name, policy_header, policy_class)
# ---------------------------------------------------------------------------
FAMILIES_PELT <- list(
  c("mean",        "mean",        "MeanFamily"),
  c("mgaussian",   "mgaussian",   "MgaussianFamily"),
  c("variance",    "variance",    "VarianceFamily"),
  c("meanvariance","meanvariance","MeanvarianceFamily"),
  c("garch",       "garch",       "GarchFamily")
)

# p == 1 specialisations: eliminate dimension loops / arma temporaries.
FAMILIES_PELT_1D <- list(
  c("mean_1d",        "mean",        "MeanFamily"),
  c("variance_1d",    "variance",    "VarianceFamily"),
  c("meanvariance_1d","meanvariance","MeanvarianceFamily")
)

FAMILIES_SEGD <- list(
  c("arma",    "arma",    "ArmaFamily"),
  c("binomial","binomial","BinomialFamily"),
  c("custom",  "custom",  "CustomFamily"),
  c("gaussian","gaussian","GaussianFamily"),
  c("lasso",   "lasso",   "LassoFamily"),
  c("ma",      "ma",      "MaFamily"),
  c("poisson", "poisson", "PoissonFamily")
)

# ---------------------------------------------------------------------------
# Variant axes  (cpp_value, suffix)
# ---------------------------------------------------------------------------
COST_ADJUSTMENTS  <- list(c("kBIC",  ""),      c("kMBIC", "_mbic"), c("kMDL", "_mdl"))
PROGRESS_VARIANTS <- list(c("false", ""),      c("true",  "_prog"))
VANILLA_VARIANTS  <- list(c("false", ""),      c("true",  "_van"))
LINE_SEARCH_VARIANTS <- list(c("false", ""),   c("true",  "_ls"))

# ---------------------------------------------------------------------------
# Generate PELT families: 5 × 3 cost × 2 progress = 30 files
# ---------------------------------------------------------------------------
for (fam in FAMILIES_PELT) {
  for (cost in COST_ADJUSTMENTS) {
    for (prog in PROGRESS_VARIANTS) {
      func_suffix <- paste0(cost[[2]], prog[[2]])
      outfile <- file.path(out_dir, paste0("fastcpd_", fam[[1]], func_suffix, ".cc"))
      generate(pelt_tmpl, list(
        "{FAMILY_NAME}"      = fam[[1]],
        "{POLICY_HEADER}"    = fam[[2]],
        "{POLICY_CLASS}"     = fam[[3]],
        "{K_R_PROGRESS}"     = prog[[1]],
        "{K_VANILLA_ONLY}"   = "false",
        "{K_COST_ADJUSTMENT}"= cost[[1]],
        "{K_LINE_SEARCH}"    = "false",
        "{K_N_DIMS}"         = "-1",
        "{FUNC_SUFFIX}"      = func_suffix
      ), outfile)
    }
  }
}

# ---------------------------------------------------------------------------
# Generate PELT p==1 specialisations: 3 × 3 cost × 2 progress = 18 files
# ---------------------------------------------------------------------------
for (fam in FAMILIES_PELT_1D) {
  for (cost in COST_ADJUSTMENTS) {
    for (prog in PROGRESS_VARIANTS) {
      func_suffix <- paste0(cost[[2]], prog[[2]])
      outfile <- file.path(out_dir, paste0("fastcpd_", fam[[1]], func_suffix, ".cc"))
      generate(pelt_tmpl, list(
        "{FAMILY_NAME}"      = fam[[1]],
        "{POLICY_HEADER}"    = fam[[2]],
        "{POLICY_CLASS}"     = fam[[3]],
        "{K_R_PROGRESS}"     = prog[[1]],
        "{K_VANILLA_ONLY}"   = "false",
        "{K_COST_ADJUSTMENT}"= cost[[1]],
        "{K_LINE_SEARCH}"    = "false",
        "{K_N_DIMS}"         = "1",
        "{FUNC_SUFFIX}"      = func_suffix
      ), outfile)
    }
  }
}

# ---------------------------------------------------------------------------
# Generate SEGD families: 7 × 3 cost × 2 van × 2 ls × 2 progress = 168 files
# ---------------------------------------------------------------------------
for (fam in FAMILIES_SEGD) {
  for (cost in COST_ADJUSTMENTS) {
    for (van in VANILLA_VARIANTS) {
      for (ls in LINE_SEARCH_VARIANTS) {
        for (prog in PROGRESS_VARIANTS) {
          func_suffix <- paste0(cost[[2]], van[[2]], ls[[2]], prog[[2]])
          outfile <- file.path(out_dir, paste0("fastcpd_", fam[[1]], func_suffix, ".cc"))
          generate(segd_tmpl, list(
            "{FAMILY_NAME}"      = fam[[1]],
            "{POLICY_HEADER}"    = fam[[2]],
            "{POLICY_CLASS}"     = fam[[3]],
            "{K_R_PROGRESS}"     = prog[[1]],
            "{K_VANILLA_ONLY}"   = van[[1]],
            "{K_COST_ADJUSTMENT}"= cost[[1]],
            "{K_LINE_SEARCH}"    = ls[[1]],
            "{K_N_DIMS}"         = "-1",
            "{FUNC_SUFFIX}"      = func_suffix
          ), outfile)
        }
      }
    }
  }
}

cat(sprintf("Generated 216 .cc files in %s\n", out_dir))
