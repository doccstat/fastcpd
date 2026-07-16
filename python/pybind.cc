#define PY_SSIZE_T_CLEAN

#include <fastcpd/fastcpd.h>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>

namespace py = pybind11;

namespace {

using DoubleArray =
    py::array_t<double, py::array::c_style | py::array::forcecast>;

arma::colvec to_colvec(DoubleArray const& values, char const* name) {
  py::buffer_info const buffer = values.request();
  if (buffer.ndim != 1) {
    throw std::invalid_argument(std::string("fastcpd: ") + name +
                                " must be one-dimensional");
  }
  arma::colvec result(static_cast<arma::uword>(buffer.shape[0]));
  if (result.n_elem > 0) {
    std::memcpy(result.memptr(), buffer.ptr,
                result.n_elem * sizeof(double));
  }
  return result;
}

arma::mat to_matrix(DoubleArray const& values, char const* name) {
  py::buffer_info const buffer = values.request();
  if (buffer.ndim != 2) {
    throw std::invalid_argument(std::string("fastcpd: ") + name +
                                " must be two-dimensional");
  }
  arma::uword const rows = static_cast<arma::uword>(buffer.shape[0]);
  arma::uword const columns = static_cast<arma::uword>(buffer.shape[1]);
  double const* const input = static_cast<double const*>(buffer.ptr);
  arma::mat result(rows, columns);
  for (arma::uword row = 0; row < rows; ++row) {
    for (arma::uword column = 0; column < columns; ++column) {
      result(row, column) = input[row * columns + column];
    }
  }
  return result;
}

py::array_t<std::int64_t> to_index_array(arma::colvec const& values) {
  py::array_t<std::int64_t> result(values.n_elem);
  std::int64_t* const output = result.mutable_data();
  for (arma::uword index = 0; index < values.n_elem; ++index) {
    output[index] = static_cast<std::int64_t>(values(index));
  }
  return result;
}

py::array_t<double> to_array(arma::colvec const& values) {
  py::array_t<double> result(values.n_elem);
  if (values.n_elem > 0) {
    std::memcpy(result.mutable_data(), values.memptr(),
                values.n_elem * sizeof(double));
  }
  return result;
}

py::array_t<double> to_array(arma::mat const& values) {
  py::array_t<double> result(
      {static_cast<py::ssize_t>(values.n_rows),
       static_cast<py::ssize_t>(values.n_cols)});
  double* const output = result.mutable_data();
  for (arma::uword row = 0; row < values.n_rows; ++row) {
    for (arma::uword column = 0; column < values.n_cols; ++column) {
      output[row * values.n_cols + column] = values(row, column);
    }
  }
  return result;
}

py::dict fastcpd_impl(
    py::object const& beta,
    std::string const& cost_adjustment,
    bool cp_only,
    DoubleArray const& data,
    double epsilon,
    std::string const& family,
    DoubleArray const& line_search,
    DoubleArray const& lower,
    double momentum_coef,
    DoubleArray const& order,
    int p,
    unsigned int p_response,
    double pruning_coef,
    unsigned int segment_count,
    double trim,
    DoubleArray const& upper,
    double vanilla_percentage,
    DoubleArray const& variance_estimate,
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
  options.line_search = to_colvec(line_search, "line_search");
  options.lower = to_colvec(lower, "lower");
  options.upper = to_colvec(upper, "upper");
  options.momentum_coef = momentum_coef;
  options.order = to_colvec(order, "order");
  options.p = p;
  options.p_response = p_response;
  if (!std::isnan(pruning_coef)) options.pruning_coef = pruning_coef;
  options.segment_count = static_cast<int>(segment_count);
  options.trim = trim;
  options.vanilla_percentage = vanilla_percentage;
  options.variance_estimate =
      to_matrix(variance_estimate, "variance_estimate");
  options.warm_start = warm_start;
  options.show_progress = show_progress;

  arma::mat data_matrix = to_matrix(data, "data");
  fastcpd::Result result;
  {
    py::gil_scoped_release release;
    result = fastcpd::detect(data_matrix, std::move(options));
  }

  py::dict output;
  output["cp_set"] = to_index_array(result.change_points);
  output["raw_cp_set"] = to_index_array(result.raw_change_points);
  output["cost_values"] = to_array(result.cost_values);
  output["residuals"] = to_array(result.residuals);
  output["thetas"] = to_array(result.thetas);
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
