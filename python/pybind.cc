// NO_RCPP must be defined before any fastcpd header so that
// fastcpd_template.h uses <armadillo> instead of <RcppArmadillo.h>.
#define NO_RCPP

#define PY_SSIZE_T_CLEAN
#define CONFIG_64
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <Python.h>

#include "fastcpd_py.h"

#include <cmath>

namespace py = pybind11;

// ---------------------------------------------------------------------------
// Internal helpers: convert Armadillo types to Python-friendly containers.
// ---------------------------------------------------------------------------

static std::vector<double> colvec_to_vec(arma::colvec const& v) {
  return std::vector<double>(v.begin(), v.end());
}

// Convert arma::mat to row-major nested vector (each inner vector is one row).
static std::vector<std::vector<double>> mat_to_rowlist(arma::mat const& m) {
  std::vector<std::vector<double>> out(m.n_rows,
                                       std::vector<double>(m.n_cols));
  for (arma::uword i = 0; i < m.n_rows; ++i)
    for (arma::uword j = 0; j < m.n_cols; ++j)
      out[i][j] = m(i, j);
  return out;
}

// ---------------------------------------------------------------------------
// Main Python-facing entry point.
//
// beta_obj  – Python str ("BIC", "MBIC", "MDL") or float.
// p         – number of model parameters; 0 triggers auto-inference in C++.
// pruning_coef – base pruning coefficient; NaN triggers auto-computation
//               (matching the default-not-set behaviour in R).
//
// Returns a Python dict with keys:
//   "cp_set"     – list of change-point indices (1-based)
//   "raw_cp_set" – raw change-point indices before boundary trimming
// and, when cp_only is False:
//   "cost_values" – list of segment cost values
//   "residuals"   – nested list (n_obs × n_response)
//   "thetas"      – nested list (n_params × n_segments), i.e. row per param
// ---------------------------------------------------------------------------

py::dict fastcpd_impl(
    py::object const& beta_obj,
    std::string const& cost_adjustment, bool const cp_only,
    std::vector<std::vector<double>> const& data, double const epsilon,
    std::string const& family, std::vector<double> const& line_search,
    std::vector<double> const& lower, double const momentum_coef,
    std::vector<double> const& order, int const p,
    unsigned int const p_response, double const pruning_coef,
    unsigned int const segment_count, double const trim,
    std::vector<double> const& upper, double const vanilla_percentage,
    std::vector<std::vector<double>> const& variance_estimate,
    bool const warm_start) {

  // -- Convert data to arma::mat (n_rows × n_cols). -------------------------
  arma::mat data_(data.size(), data.empty() ? 0 : data[0].size());
  for (size_t i = 0; i < data.size(); ++i)
    for (size_t j = 0; j < data[i].size(); ++j)
      data_(i, j) = data[i][j];

  // -- Convert 1-D vectors to arma::colvec. ---------------------------------
  auto to_colvec = [](std::vector<double> const& v) {
    arma::colvec c(v.size());
    for (size_t i = 0; i < v.size(); ++i) c(i) = v[i];
    return c;
  };
  arma::colvec line_search_ = to_colvec(
      line_search.empty() ? std::vector<double>{1.0} : line_search);
  arma::colvec lower_  = to_colvec(lower);
  arma::colvec order_  = to_colvec(order);
  arma::colvec upper_  = to_colvec(upper);

  // -- Convert variance_estimate to arma::mat. ------------------------------
  arma::mat variance_estimate_;
  if (variance_estimate.empty() || variance_estimate[0].empty()) {
    variance_estimate_ = arma::eye(1, 1);
  } else {
    variance_estimate_.set_size(variance_estimate.size(),
                                variance_estimate[0].size());
    for (size_t i = 0; i < variance_estimate.size(); ++i)
      for (size_t j = 0; j < variance_estimate[i].size(); ++j)
        variance_estimate_(i, j) = variance_estimate[i][j];
  }

  // -- Resolve p, beta, pruning_coef in C++ ---------------------------------
  int const n_rows = static_cast<int>(data_.n_rows);
  int const n_cols = static_cast<int>(data_.n_cols);

  int const resolved_p = fastcpd::fastcpd_compute_p(
      family, n_cols, static_cast<int>(p_response), p);

  // beta_obj is either a Python str ("MBIC" etc.) or a numeric float.
  double beta_val;
  std::string beta_str;
  if (py::isinstance<py::str>(beta_obj)) {
    beta_str = beta_obj.cast<std::string>();
    beta_val = 0.0;  // unused; fastcpd_compute_beta ignores when str non-empty
  } else {
    beta_val = beta_obj.cast<double>();
  }
  double const resolved_beta = fastcpd::fastcpd_compute_beta(
      beta_val, beta_str, n_rows, resolved_p);

  // NaN pruning_coef signals "not explicitly set" – apply auto-adjustment.
  bool const pruning_explicitly_set = !std::isnan(pruning_coef);
  double const base_pruning = pruning_explicitly_set ? pruning_coef : 0.0;
  double const resolved_pruning = fastcpd::fastcpd_compute_pruning_coef(
      base_pruning, pruning_explicitly_set, cost_adjustment, family,
      resolved_p);

  // -- Run the C++ detection. -----------------------------------------------
  auto result = fastcpd::fastcpd_py_dispatch(
      resolved_beta, cost_adjustment, cp_only, data_, epsilon, family,
      line_search_, lower_, momentum_coef, order_,
      resolved_p, p_response, resolved_pruning,
      static_cast<int>(segment_count), trim, upper_,
      vanilla_percentage, variance_estimate_, warm_start);

  // -- Build the Python dict. -----------------------------------------------
  py::dict out;
  out["cp_set"]     = colvec_to_vec(std::get<1>(result));
  out["raw_cp_set"] = colvec_to_vec(std::get<0>(result));
  if (!cp_only) {
    out["cost_values"] = colvec_to_vec(std::get<2>(result));
    out["residuals"]   = mat_to_rowlist(std::get<3>(result));
    out["thetas"]      = mat_to_rowlist(std::get<4>(result));
  }
  return out;
}

PYBIND11_MODULE(interface, module) {
  module.doc() = "fastcpd C++/Python interface";
  module.def("fastcpd_impl", &fastcpd_impl,
             "Fast change-point detection");
}
