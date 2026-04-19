#pragma once
// Python-facing dispatch header (NO_RCPP mode).
// Included by python/pybind.cc; must not include any Rcpp headers.

#include <armadillo>
#include <string>
#include <tuple>
#include <vector>

namespace fastcpd {

// ---------------------------------------------------------------------------
// Helper free functions – called by pybind.cc before dispatching.
// These mirror the parameter-normalisation logic in R/fastcpd.R and
// R/utilities.R so that neither the R nor Python wrapper needs to duplicate
// the formulas.
// ---------------------------------------------------------------------------

/// Infer the number of model parameters from the family and data dimensions.
/// Returns p_explicit when p_explicit > 0; auto-infers from family otherwise.
///
/// Mirrors the p-computation blocks in R/fastcpd.R lines 303-375 and
/// Python segmentation.py lines 207-231 for the families supported by the
/// Python binding.
int fastcpd_compute_p(std::string const& family, int n_cols,
                      int p_response, int p_explicit = 0);

/// Convert a string beta criterion ("BIC", "MBIC", "MDL") to a numeric
/// penalty value.  When beta_str is empty the supplied beta_val is returned
/// directly.
///
/// Mirrors R/fastcpd.R lines 401-419 (without the gaussian σ scaling, which
/// is R-only) and Python segmentation.py lines 241-249.
double fastcpd_compute_beta(double beta_val, std::string const& beta_str,
                            int n, int p);

/// Compute the final pruning coefficient.
///
/// When pruning_coef_explicitly_set is false:
///   – sets the base to –∞ for families that must never prune (mgaussian,
///     lasso), matching R's get_pruning_coef() in utilities.R lines 127-136;
///   – adds p·log(2) for MBIC or p·log₂(2) for MDL.
/// When pruning_coef_explicitly_set is true the supplied value is returned
/// unchanged.
double fastcpd_compute_pruning_coef(double pruning_coef,
                                    bool pruning_coef_explicitly_set,
                                    std::string const& cost_adjustment,
                                    std::string const& family, int p);

// ---------------------------------------------------------------------------
// Main dispatch entry point.
// ---------------------------------------------------------------------------

/// Dispatch change-point detection to the correct Fastcpd<> template.
/// Returns (raw_cp_set, cp_set, cost_values, residual, thetas).
///
/// Supported families: mean, variance, meanvariance, mgaussian, lasso.
/// (Families requiring Rcpp — arma, ma, garch, gaussian/lm, binomial, poisson
///  — are not yet supported in the Python binding.)
std::tuple<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat>
fastcpd_py_dispatch(
    double beta, std::string const& cost_adjustment, bool cp_only,
    arma::mat const& data, double epsilon, std::string const& family,
    arma::colvec const& line_search, arma::colvec const& lower,
    double momentum_coef, arma::colvec const& order, int p,
    unsigned int p_response, double pruning_coef, int segment_count,
    double trim, arma::colvec const& upper, double vanilla_percentage,
    arma::mat const& variance_estimate, bool warm_start);

}  // namespace fastcpd
