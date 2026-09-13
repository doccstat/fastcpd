#define PY_SSIZE_T_CLEAN

#include <fastcpd/fastcpd.h>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace py = pybind11;

namespace {

using DoubleArray =
    py::array_t<double, py::array::c_style | py::array::forcecast>;

arma::colvec to_colvec(DoubleArray const& values, char const* name);
arma::mat to_matrix(DoubleArray const& values, char const* name);

// Python callbacks are invoked while the detector's outer GIL-release scope
// is active. Keep the Python function in a shared holder whose final
// destruction reacquires the GIL; native solver copies of std::function may
// otherwise release the last Python reference on a GIL-free thread.
using Callback = std::shared_ptr<py::function>;

Callback make_callback(py::object const& object, char const* name) {
  if (object.is_none()) return {};
  if (!PyCallable_Check(object.ptr())) {
    throw py::type_error(std::string("fastcpd: ") + name +
                         " must be callable");
  }
  return Callback(
      new py::function(object), [](py::function* function) {
        py::gil_scoped_acquire acquire;
        delete function;
      });
}

double callback_scalar(py::object const& value, char const* name) {
  double result;
  try {
    result = value.cast<double>();
  } catch (py::cast_error const& error) {
    throw py::type_error(std::string("fastcpd: ") + name +
                         " must return a numeric scalar");
  }
  // Infinite costs are useful for custom models that disallow short or
  // otherwise invalid segments. NaN has no meaningful ordering in PELT/SEN.
  if (std::isnan(result)) {
    throw std::invalid_argument(std::string("fastcpd: ") + name +
                                " must not return NaN");
  }
  return result;
}

arma::colvec callback_colvec(py::object const& value, char const* name,
                             arma::uword expected_size) {
  DoubleArray array = DoubleArray::ensure(value);
  if (!array) {
    throw py::type_error(std::string("fastcpd: ") + name +
                         " must return a one-dimensional numeric array");
  }
  arma::colvec result = to_colvec(array, name);
  if (result.n_elem != expected_size) {
    throw std::invalid_argument(std::string("fastcpd: ") + name +
                                " returned " +
                                std::to_string(result.n_elem) +
                                " values; expected " +
                                std::to_string(expected_size));
  }
  if (!result.is_finite()) {
    throw std::invalid_argument(std::string("fastcpd: ") + name +
                                " must return finite values");
  }
  return result;
}

arma::mat callback_matrix(py::object const& value, char const* name,
                          arma::uword expected_size) {
  DoubleArray array = DoubleArray::ensure(value);
  if (!array) {
    throw py::type_error(std::string("fastcpd: ") + name +
                         " must return a two-dimensional numeric array");
  }
  arma::mat result = to_matrix(array, name);
  if (result.n_rows != expected_size || result.n_cols != expected_size) {
    throw std::invalid_argument(std::string("fastcpd: ") + name +
                                " returned a " +
                                std::to_string(result.n_rows) + "x" +
                                std::to_string(result.n_cols) +
                                " matrix; expected " +
                                std::to_string(expected_size) + "x" +
                                std::to_string(expected_size));
  }
  if (!result.is_finite()) {
    throw std::invalid_argument(std::string("fastcpd: ") + name +
                                " must return finite values");
  }
  return result;
}

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
    bool show_progress,
    py::object cost_pelt,
    py::object cost_sen,
    py::object cost_gradient,
    py::object cost_hessian) {
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

  Callback const pelt_callback = make_callback(cost_pelt, "cost_pelt");
  Callback const sen_callback = make_callback(cost_sen, "cost_sen");
  Callback const gradient_callback =
      make_callback(cost_gradient, "cost_gradient");
  Callback const hessian_callback =
      make_callback(cost_hessian, "cost_hessian");
  if (pelt_callback) {
    options.cost_pelt = [pelt_callback](arma::mat const& segment) {
      py::gil_scoped_acquire acquire;
      py::object const value = (*pelt_callback)(to_array(segment));
      return callback_scalar(value, "cost_pelt");
    };
  }
  if (sen_callback) {
    options.cost_sen = [sen_callback](arma::mat const& segment,
                                      arma::colvec const& theta) {
      py::gil_scoped_acquire acquire;
      py::object const value =
          (*sen_callback)(to_array(segment), to_array(theta));
      return callback_scalar(value, "cost_sen");
    };
  }
  if (gradient_callback) {
    options.cost_gradient =
        [gradient_callback](arma::mat const& segment,
                            arma::colvec const& theta) {
          py::gil_scoped_acquire acquire;
          py::object const value =
              (*gradient_callback)(to_array(segment), to_array(theta));
          return callback_colvec(value, "cost_gradient", theta.n_elem);
        };
  }
  if (hessian_callback) {
    options.cost_hessian =
        [hessian_callback](arma::mat const& segment,
                           arma::colvec const& theta) {
          py::gil_scoped_acquire acquire;
          py::object const value =
              (*hessian_callback)(to_array(segment), to_array(theta));
          return callback_matrix(value, "cost_hessian", theta.n_elem);
        };
  }

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
      py::arg("show_progress") = false,
      py::arg("cost_pelt") = py::none(), py::arg("cost_sen") = py::none(),
      py::arg("cost_gradient") = py::none(),
      py::arg("cost_hessian") = py::none());
}
