#define PY_SSIZE_T_CLEAN

#include <fastcpd/fastcpd.h>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace py = pybind11;

namespace {

arma::colvec to_colvec(std::vector<double> const& values) {
  arma::colvec result(values.size());
  for (std::size_t index = 0; index < values.size(); ++index) {
    result(index) = values[index];
  }
  return result;
}

arma::mat to_matrix(std::vector<std::vector<double>> const& values,
                    char const* name) {
  if (values.empty()) return arma::mat();
  std::size_t const columns = values.front().size();
  if (columns == 0) return arma::mat(values.size(), 0);
  for (auto const& row : values) {
    if (row.size() != columns) {
      throw std::invalid_argument(std::string("fastcpd: ") + name +
                                  " must be rectangular");
    }
  }

  arma::mat result(values.size(), columns);
  for (std::size_t row = 0; row < values.size(); ++row) {
    for (std::size_t column = 0; column < columns; ++column) {
      result(row, column) = values[row][column];
    }
  }
  return result;
}

std::vector<double> to_vector(arma::colvec const& values) {
  return std::vector<double>(values.begin(), values.end());
}

std::vector<std::vector<double>> to_rows(arma::mat const& values) {
  std::vector<std::vector<double>> result(
      values.n_rows, std::vector<double>(values.n_cols));
  for (arma::uword row = 0; row < values.n_rows; ++row) {
    for (arma::uword column = 0; column < values.n_cols; ++column) {
      result[row][column] = values(row, column);
    }
  }
  return result;
}

py::dict fastcpd_impl(
    py::object const& beta,
    std::string const& cost_adjustment,
    bool cp_only,
    std::vector<std::vector<double>> const& data,
    double epsilon,
    std::string const& family,
    std::vector<double> const& line_search,
    std::vector<double> const& lower,
    double momentum_coef,
    std::vector<double> const& order,
    int p,
    unsigned int p_response,
    double pruning_coef,
    unsigned int segment_count,
    double trim,
    std::vector<double> const& upper,
    double vanilla_percentage,
    std::vector<std::vector<double>> const& variance_estimate,
    bool warm_start,
    bool show_progress) {
  fastcpd::Options options;
  options.family = family;
  if (py::isinstance<py::str>(beta)) {
    options.beta_criterion = beta.cast<std::string>();
  } else {
    options.beta = beta.cast<double>();
  }
  options.cost_adjustment = cost_adjustment;
  options.cp_only = cp_only;
  options.epsilon = epsilon;
  options.line_search = to_colvec(line_search);
  options.lower = to_colvec(lower);
  options.upper = to_colvec(upper);
  options.momentum_coef = momentum_coef;
  options.order = to_colvec(order);
  options.p = p;
  options.p_response = p_response;
  if (!std::isnan(pruning_coef)) options.pruning_coef = pruning_coef;
  options.segment_count = static_cast<int>(segment_count);
  options.trim = trim;
  options.vanilla_percentage = vanilla_percentage;
  options.variance_estimate = to_matrix(variance_estimate, "variance_estimate");
  options.warm_start = warm_start;
  options.show_progress = show_progress;

  fastcpd::Result const result =
      fastcpd::detect(to_matrix(data, "data"), std::move(options));

  py::dict output;
  output["cp_set"] = to_vector(result.change_points);
  output["raw_cp_set"] = to_vector(result.raw_change_points);
  if (!cp_only) {
    output["cost_values"] = to_vector(result.cost_values);
    output["residuals"] = to_rows(result.residuals);
    output["thetas"] = to_rows(result.thetas);
  }
  return output;
}

}  // namespace

PYBIND11_MODULE(interface, module) {
  module.doc() =
      "Python bindings for the shared fastcpd standalone C++ implementation";
  module.def(
      "fastcpd_impl", &fastcpd_impl, "Fast change-point detection",
      py::arg("beta"), py::arg("cost_adjustment"), py::arg("cp_only"),
      py::arg("data"), py::arg("epsilon"), py::arg("family"),
      py::arg("line_search"), py::arg("lower"), py::arg("momentum_coef"),
      py::arg("order"), py::arg("p"), py::arg("p_response"),
      py::arg("pruning_coef"), py::arg("segment_count"), py::arg("trim"),
      py::arg("upper"), py::arg("vanilla_percentage"),
      py::arg("variance_estimate"), py::arg("warm_start"),
      py::arg("show_progress") = false);
}
