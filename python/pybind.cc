#include "src/fastcpd.h"

#define PY_SSIZE_T_CLEAN
#define CONFIG_64
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <Python.h>

std::vector<double> fastcpd_impl(
    double const beta, std::string const& cost_adjustment, bool const cp_only,
    std::vector<std::vector<double>> const& data, double const epsilon,
    std::string const& family, std::vector<double> const& line_search,
    std::vector<double> const& lower, double const momentum_coef,
    std::vector<double> const& order, unsigned int const p,
    unsigned int const p_response, double const pruning_coef,
    unsigned int const segment_count, double const trim,
    std::vector<double> const& upper, double const vanilla_percentage,
    std::vector<std::vector<double>> const& variance_estimate,
    bool const warm_start) {
  // Convert data to arma::mat
  arma::mat data_(data.size(), data[0].size());
  for (size_t i = 0; i < data.size(); ++i) {
    for (size_t j = 0; j < data[i].size(); ++j) {
      data_(i, j) = data[i][j];
    }
  }

  // Convert line_search to arma::colvec
  arma::colvec line_search_(line_search.size());
  for (size_t i = 0; i < line_search.size(); ++i) {
    line_search_(i) = line_search[i];
  }

  // Convert lower to arma::colvec
  arma::colvec lower_(lower.size());
  for (size_t i = 0; i < lower.size(); ++i) {
    lower_(i) = lower[i];
  }

  // Convert order to arma::colvec
  arma::colvec order_(order.size());
  for (size_t i = 0; i < order.size(); ++i) {
    order_(i) = order[i];
  }

  // Convert upper to arma::colvec
  arma::colvec upper_(upper.size());
  for (size_t i = 0; i < upper.size(); ++i) {
    upper_(i) = upper[i];
  }

  // Convert variance_estimate to arma::mat
  arma::mat variance_estimate_(variance_estimate.size(),
                               variance_estimate[0].size());
  for (size_t i = 0; i < variance_estimate.size(); ++i) {
    for (size_t j = 0; j < variance_estimate[i].size(); ++j) {
      variance_estimate_(i, j) = variance_estimate[i][j];
    }
  }

  fastcpd::classes::Fastcpd fastcpd_instance(
      /* beta */ beta,
      /* cost_pelt */ nullptr,
      /* cost_sen */ nullptr,
      /* cost_adjustment */ cost_adjustment,
      /* cost_gradient */ nullptr,
      /* cost_hessian */ nullptr,
      /* cp_only */ cp_only,
      /* data */ data_,
      /* epsilon */ epsilon,
      /* family */ family,
      /* multiple_epochs_function */ nullptr,
      /* line_search */ line_search_,
      /* lower */ lower_,
      /* momentum_coef */ momentum_coef,
      /* order */ order_,
      /* p */ p,
      /* p_response */ p_response,
      /* pruning_coef */ pruning_coef,
      /* r_progress */ false,
      /* segment_count */ segment_count,
      /* trim */ trim,
      /* upper */ upper_,
      /* vanilla_percentage */ vanilla_percentage,
      /* variance_estimate */ variance_estimate_,
      /* warm_start */ warm_start);

  std::tuple<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat>
      result = fastcpd_instance.Run();

  arma::colvec change_points = std::get<1>(result);

  return std::vector<double>(change_points.begin(), change_points.end());
}

PYBIND11_MODULE(interface, module) {
  module.doc() = "fastcpd C++/Python interface";
  module.def("fastcpd_impl", &fastcpd_impl,
             "A function that computes the fast change point detection");
}
