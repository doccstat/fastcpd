#ifndef NO_RCPP
#define NO_RCPP
#endif

#include <fastcpd/fastcpd.h>

#include "fastcpd_template.h"
#include "families/arima.h"
#include "families/arma.h"
#include "families/binomial.h"
#include "families/custom.h"
#include "families/exponential.h"
#include "families/garch.h"
#include "families/gaussian.h"
#include "families/lasso.h"
#include "families/ma.h"
#include "families/mean.h"
#include "families/meanvariance.h"
#include "families/mgaussian.h"
#include "families/poisson.h"
#include "families/quantile.h"
#include "families/variance.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace fastcpd {
namespace detail {

using RunResult =
    std::tuple<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat>;

std::string lower_ascii(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char ch) {
                   return static_cast<char>(std::tolower(ch));
                 });
  return value;
}

bool is_pelt_family(std::string const& family) {
  return family == "mean" || family == "variance" ||
         family == "meanvariance" || family == "exponential" ||
         family == "mgaussian" || family == "garch" || family == "arima";
}

void validate_common_options(Options const& options) {
  auto require_finite = [](double value, char const* name) {
    if (!std::isfinite(value)) {
      throw std::invalid_argument(std::string("fastcpd: ") + name +
                                  " must be finite");
    }
  };

  if (!std::isfinite(options.trim) || options.trim < 0.0 ||
      options.trim > 1.0) {
    throw std::invalid_argument("fastcpd: trim must be in [0, 1]");
  }
  if (!std::isfinite(options.vanilla_percentage) ||
      options.vanilla_percentage < 0.0 ||
      options.vanilla_percentage > 1.0) {
    throw std::invalid_argument(
        "fastcpd: vanilla_percentage must be in [0, 1]");
  }
  if (options.segment_count <= 0) {
    throw std::invalid_argument(
        "fastcpd: segment_count must be a positive integer");
  }
  if (!std::isfinite(options.epsilon) || options.epsilon <= 0.0) {
    throw std::invalid_argument("fastcpd: epsilon must be finite and positive");
  }
  require_finite(options.momentum_coef, "momentum_coef");
  if (options.p < 0) {
    throw std::invalid_argument("fastcpd: p must be non-negative");
  }
  if (options.beta.has_value()) require_finite(*options.beta, "beta");
  if (options.pruning_coef.has_value()) {
    require_finite(*options.pruning_coef, "pruning_coef");
  }

  for (double const value : options.line_search) {
    if (!std::isfinite(value) || value <= 0.0) {
      throw std::invalid_argument(
          "fastcpd: line_search values must be finite and positive");
    }
  }
  for (double const value : options.lower) {
    if (std::isnan(value)) {
      throw std::invalid_argument("fastcpd: lower must not contain NaN");
    }
  }
  for (double const value : options.upper) {
    if (std::isnan(value)) {
      throw std::invalid_argument("fastcpd: upper must not contain NaN");
    }
  }
  if (options.lower.n_elem > 0 && options.upper.n_elem > 0 &&
      options.lower.n_elem != options.upper.n_elem) {
    throw std::invalid_argument(
        "fastcpd: lower and upper must have the same length");
  }
  if (options.lower.n_elem > 0 && options.upper.n_elem > 0 &&
      arma::any(options.lower > options.upper)) {
    throw std::invalid_argument(
        "fastcpd: lower values must not exceed upper values");
  }
}

void validate_cost_adjustment(std::string const& value) {
  if (value != "BIC" && value != "MBIC" && value != "MDL") {
    throw std::invalid_argument(
        "fastcpd: cost_adjustment must be BIC, MBIC, or MDL");
  }
}

std::string normalize_family(std::string family, arma::colvec const& order) {
  family = lower_ascii(std::move(family));
  if (family == "lm") return "gaussian";
  if (family == "var") return "mgaussian";
  if (family == "arma" && order.n_elem >= 1 && order(0) == 0.0) return "ma";
  return family;
}

bool is_nonnegative_integer(double value) {
  return std::isfinite(value) && value >= 0.0 && value == std::floor(value);
}

void validate_integer_order(arma::colvec const& order,
                            arma::uword expected_length,
                            char const* family) {
  if (order.n_elem != expected_length) {
    throw std::invalid_argument(std::string("fastcpd: ") + family +
                                " order has the wrong length");
  }
  for (double const value : order) {
    if (!is_nonnegative_integer(value)) {
      throw std::invalid_argument(std::string("fastcpd: ") + family +
                                  " order values must be non-negative integers");
    }
  }
}

void validate_family_input(std::string const& family, arma::mat const& data,
                           arma::colvec const& order,
                           unsigned int p_response) {
  if (!data.is_finite()) {
    throw std::invalid_argument(
        "fastcpd: data must contain only finite values");
  }
  bool const regression_family =
      family == "gaussian" || family == "lasso" || family == "binomial" ||
      family == "poisson" || family == "quantile";
  if (regression_family && data.n_cols < 2) {
    throw std::invalid_argument(
        "fastcpd: regression data must include predictors");
  }
  if (family == "mgaussian") {
    if (p_response == 0 || p_response >= data.n_cols) {
      throw std::invalid_argument(
          "fastcpd: mgaussian requires p_response leading response columns "
          "and at least one predictor column");
    }
  } else if (p_response > data.n_cols) {
    throw std::invalid_argument(
        "fastcpd: p_response must not exceed data.n_cols");
  }

  if (family == "arma") {
    validate_integer_order(order, 2, "ARMA");
    if (order(0) == 0.0) {
      throw std::invalid_argument(
          "fastcpd: ARMA order must have a positive AR component after "
          "normalization");
    }
  } else if (family == "ma") {
    validate_integer_order(order, 2, "MA");
    if (order(0) != 0.0 || order(1) <= 0.0) {
      throw std::invalid_argument(
          "fastcpd: MA order must be (0, q), with q positive");
    }
  } else if (family == "garch") {
    validate_integer_order(order, 2, "GARCH");
    if (order(0) == 0.0 && order(1) == 0.0) {
      throw std::invalid_argument(
          "fastcpd: GARCH order must contain a non-zero value");
    }
  } else if (family == "arima") {
    validate_integer_order(order, 3, "ARIMA");
    if (arma::all(order == 0.0)) {
      throw std::invalid_argument(
          "fastcpd: ARIMA order must contain a non-zero value");
    }
    if (order(1) >= static_cast<double>(data.n_rows)) {
      throw std::invalid_argument(
          "fastcpd: ARIMA integration order must be smaller than data.n_rows");
    }
  } else if (family == "quantile") {
    if (order.n_elem != 1 || !std::isfinite(order(0)) || order(0) <= 0.0 ||
        order(0) >= 1.0) {
      throw std::invalid_argument(
          "fastcpd: quantile order must contain one value in (0, 1)");
    }
  }

  if ((family == "arma" || family == "ma" || family == "garch" ||
       family == "arima") &&
      data.n_cols != 1) {
    throw std::invalid_argument(std::string("fastcpd: ") + family +
                                " data must be univariate");
  }
}

double order_at(arma::colvec const& order, arma::uword index,
                double fallback = 0.0) {
  return index < order.n_elem ? order(index) : fallback;
}

int compute_p(std::string const& family, arma::mat const& data,
              arma::colvec const& order, unsigned int p_response,
              int p_explicit) {
  if (p_explicit > 0) return p_explicit;
  int const n_cols = static_cast<int>(data.n_cols);
  if (family == "mean" || family == "exponential") return n_cols;
  if (family == "variance") return n_cols * n_cols;
  if (family == "meanvariance") return n_cols + n_cols * n_cols;
  if (family == "mgaussian") {
    int const q = p_response > 0 ? static_cast<int>(p_response) : n_cols;
    int const predictor_count = n_cols - q;
    return predictor_count > 0 ? q * predictor_count : q;
  }
  if (family == "gaussian" || family == "lasso" ||
      family == "binomial" || family == "poisson" ||
      family == "quantile") {
    return n_cols - 1;
  }
  if (family == "garch" || family == "arma") {
    return static_cast<int>(arma::sum(order)) + 1;
  }
  if (family == "arima") {
    if (order.n_elem != 3) {
      throw std::invalid_argument(
          "fastcpd: ARIMA order must contain p, d, and q");
    }
    bool any_nonzero = false;
    for (arma::uword index = 0; index < order.n_elem; ++index) {
      double const value = order(index);
      if (!std::isfinite(value) || value < 0.0 || value != std::floor(value)) {
        throw std::invalid_argument(
            "fastcpd: ARIMA order values must be non-negative integers");
      }
      any_nonzero = any_nonzero || value != 0.0;
    }
    if (!any_nonzero) {
      throw std::invalid_argument(
          "fastcpd: ARIMA order must contain a non-zero value");
    }
    return static_cast<int>(order(0) + order(2)) + 1;
  }
  if (family == "ma") return static_cast<int>(order_at(order, 1)) + 1;
  if (family == "custom") return std::max(1, n_cols - 1);
  throw std::invalid_argument("fastcpd: unsupported family '" + family + "'");
}

double compute_beta(std::optional<double> beta, std::string const& criterion,
                    int n, int p) {
  if (beta.has_value()) return *beta;
  if (criterion == "BIC") {
    return (p + 1) * std::log(static_cast<double>(n)) / 2.0;
  }
  if (criterion == "MBIC") {
    return (p + 2) * std::log(static_cast<double>(n)) / 2.0;
  }
  if (criterion == "MDL") {
    return (p + 2) * std::log2(static_cast<double>(n)) / 2.0;
  }
  throw std::invalid_argument(
      "fastcpd: beta_criterion must be BIC, MBIC, or MDL");
}

double compute_pruning_coef(std::optional<double> pruning_coef,
                            std::string const& cost_adjustment,
                            std::string const& family, int p) {
  double value = pruning_coef.value_or(0.0);
  if (!pruning_coef.has_value() &&
      (family == "mgaussian" || family == "lasso" || family == "garch" ||
       family == "arima")) {
    value = -std::numeric_limits<double>::infinity();
  }
  if (!pruning_coef.has_value() && cost_adjustment == "MBIC") {
    value += p * std::log(2.0);
  } else if (!pruning_coef.has_value() && cost_adjustment == "MDL") {
    value += p * std::log2(2.0);
  }
  return value;
}

arma::colvec fill_or_validate_bound(arma::colvec const& value, int p,
                                    double fill, char const* name) {
  if (value.n_elem == 0) {
    arma::colvec out(static_cast<arma::uword>(p));
    out.fill(fill);
    return out;
  }
  if (value.n_elem != static_cast<arma::uword>(p)) {
    throw std::invalid_argument(std::string("fastcpd: ") + name +
                                " must have length p");
  }
  return value;
}

arma::colvec normalize_line_search(arma::colvec const& value) {
  if (value.n_elem == 0) return arma::colvec{1.0};
  return value;
}

// Match the R helper `nearest_pd_`: symmetrize through an eigendecomposition
// and lift non-positive eigenvalues only to machine epsilon.  R applies this
// projection to internally estimated multivariate-regression covariance
// matrices, but not to user-supplied matrices (see below).
arma::mat nearest_positive_definite(arma::mat const& matrix) {
  if (matrix.n_rows != matrix.n_cols) {
    throw std::invalid_argument(
        "fastcpd: variance_estimate must be a square matrix");
  }
  if (matrix.n_rows == 0) return matrix;
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  // R's symmetric eigen path consumes the lower triangle (`uplo = "L"`).
  // Use symmatl rather than Armadillo's upper-triangle symmatu so asymmetric
  // block estimates follow the same convention.
  arma::eig_sym(eigenvalues, eigenvectors, arma::symmatl(matrix));
  // `.Machine$double.eps` is the exact floor used by R's nearest_pd_.
  double const floor = std::numeric_limits<double>::epsilon();
  eigenvalues.transform([floor](double value) {
    return value > floor ? value : floor;
  });
  return eigenvectors * arma::diagmat(eigenvalues) * eigenvectors.t();
}

// Apply the common R `rcond(sigma_) < 1e-10` guard.  A well-conditioned
// user-supplied covariance is returned byte-for-byte unchanged; R replaces an
// ill-conditioned one with 1e-10 * I rather than projecting it to a nearby
// SPD matrix. R's rcond() errors for NaN/Inf matrices, so reject those rather
// than silently substituting a covariance that changes the requested fit.
arma::mat apply_rcond_fallback(arma::mat value) {
  if (value.n_rows != value.n_cols) {
    throw std::invalid_argument(
        "fastcpd: variance_estimate must be a square matrix");
  }
  if (value.n_rows == 0) return value;
  if (!value.is_finite()) {
    throw std::invalid_argument(
        "fastcpd: variance_estimate must contain only finite values");
  }

  double const reciprocal_condition = arma::rcond(value);
  if (!std::isfinite(reciprocal_condition)) {
    throw std::runtime_error(
        "fastcpd: unable to compute variance_estimate condition number");
  }
  if (reciprocal_condition < 1e-10) {
    value = 1e-10 * arma::eye<arma::mat>(value.n_rows, value.n_cols);
  }
  return value;
}

arma::mat estimate_mean_variance(arma::mat const& data) {
  if (data.n_rows < 2) return arma::eye(data.n_cols, data.n_cols);
  arma::mat const diffs =
      data.rows(1, data.n_rows - 1) - data.rows(0, data.n_rows - 2);
  // R's variance.mean() returns the raw Rice estimate.  The common rcond
  // guard (applied by resolve_variance_estimate) handles singular estimates.
  return diffs.t() * diffs /
         (2.0 * static_cast<double>(diffs.n_rows));
}

// Reproduce R's block-lagged `variance.lm()` estimator.  `data` is laid out
// as [responses, predictors], with `d` response columns.  The default block
// size in R is ncol(data) - d + 1, so adjacent blocks share all but one row.
// Each block contributes an element-wise d×d estimate; means are taken
// independently while ignoring NA values.
arma::mat estimate_linear_regression_variance(arma::mat const& data,
                                              unsigned int d,
                                              unsigned int block_size = 0,
                                              double outlier_iqr =
                                                  std::numeric_limits<
                                                      double>::infinity()) {
  if (d == 0 || d > data.n_cols) {
    throw std::invalid_argument(
        "fastcpd: p_response must be between 1 and data.n_cols");
  }
  unsigned int const predictor_count =
      static_cast<unsigned int>(data.n_cols - d);
  // R's default block size is predictor_count + 1.  With no predictors the
  // regression solve is 0×0 and variance.lm() yields no usable estimates;
  // callers handle this degenerate direct-mgaussian case below.
  if (block_size == 0) block_size = predictor_count + 1;
  if (block_size == 0 || data.n_rows <= block_size) {
    throw std::invalid_argument(
        "fastcpd: block_size must be positive and smaller than data.n_rows");
  }
  if (predictor_count == 0) {
    return arma::mat(d, d, arma::fill::value(
                              std::numeric_limits<double>::quiet_NaN()));
  }

  arma::mat sums(d, d, arma::fill::zeros);
  arma::umat counts(d, d, arma::fill::zeros);
  std::vector<double> scalar_estimates;
  arma::mat const y = data.cols(0, d - 1);
  arma::mat const x = data.cols(d, data.n_cols - 1);

  auto reciprocal_condition_ok = [](arma::mat const& value) -> bool {
    if (!value.is_finite()) return false;
    try {
      double const reciprocal_condition = arma::rcond(value);
      return std::isfinite(reciprocal_condition) &&
             reciprocal_condition >= std::numeric_limits<double>::epsilon();
    } catch (...) {
      return false;
    }
  };
  auto solve_square = [](arma::mat const& lhs, arma::mat const& rhs,
                         arma::mat* result) -> bool {
    // The status overload avoids warning/throw paths and mirrors R's
    // tryCatch around each block.  `no_approx` keeps singular blocks as
    // failures instead of silently using a least-squares approximation.
    bool const ok = arma::solve(*result, lhs, rhs, arma::solve_opts::no_approx);
    return ok && result->is_finite();
  };
  auto invert_square = [](arma::mat const& matrix, arma::mat* result) -> bool {
    bool const ok = arma::inv(*result, matrix);
    return ok && result->is_finite();
  };

  arma::uword const block_count = data.n_rows - block_size;
  for (arma::uword i = 0; i < block_count; ++i) {
    arma::mat const x_block = x.rows(i, i + block_size - 1);
    arma::mat const x_block_lagged = x.rows(i + 1, i + block_size);
    arma::mat const y_block = y.rows(i, i + block_size - 1);
    arma::mat const y_block_lagged = y.rows(i + 1, i + block_size);
    arma::mat const x_t_x = x_block.t() * x_block;
    arma::mat const x_t_x_lagged = x_block_lagged.t() * x_block_lagged;

    arma::mat block_slope;
    arma::mat block_lagged_slope;
    arma::mat x_t_x_inv;
    arma::mat x_t_x_inv_lagged;
    if (!reciprocal_condition_ok(x_t_x) ||
        !reciprocal_condition_ok(x_t_x_lagged) ||
        !solve_square(x_t_x, x_block.t() * y_block, &block_slope) ||
        !solve_square(x_t_x_lagged, x_block_lagged.t() * y_block_lagged,
                      &block_lagged_slope) ||
        !invert_square(x_t_x, &x_t_x_inv) ||
        !invert_square(x_t_x_lagged, &x_t_x_inv_lagged)) {
      // R's tryCatch stores an all-NA block when any solve fails.
      continue;
    }

    arma::mat const cross_term_x =
        x_block.rows(1, block_size - 1).t() *
        x_block_lagged.rows(0, block_size - 2);
    arma::mat const cross_term =
        x_t_x_inv * x_t_x_inv_lagged * cross_term_x;
    arma::mat const slope_delta = block_slope - block_lagged_slope;
    arma::mat const delta_numerator = slope_delta.t() * slope_delta;
    arma::mat delta_denominator(d, d, arma::fill::zeros);
    for (unsigned int j = 0; j < d; ++j) {
      for (unsigned int k = 0; k < d; ++k) {
        if (j != k) {
          arma::colvec const delta =
              block_slope.col(j) - block_lagged_slope.col(k);
          delta_denominator(j, k) = arma::dot(delta, delta);
        }
      }
    }
    // R adds this scalar to every d×d entry via vector recycling.
    double const denominator_offset =
        arma::trace(x_t_x_inv + x_t_x_inv_lagged - 2.0 * cross_term);
    delta_denominator += denominator_offset;
    arma::mat const estimate = delta_numerator / delta_denominator;

    for (unsigned int j = 0; j < d; ++j) {
      for (unsigned int k = 0; k < d; ++k) {
        double const cell = estimate(j, k);
        // `na.rm = TRUE` removes NA/NaN but retains +/-Inf.  Infinite values
        // are not expected for valid blocks, but preserving them is closer to
        // R than silently dropping all non-finite values.
        if (!std::isnan(cell)) {
          if (d == 1) {
            scalar_estimates.push_back(cell);
          } else {
            sums(j, k) += cell;
            counts(j, k) += 1u;
          }
        }
      }
    }
  }

  arma::mat result(d, d, arma::fill::zeros);
  if (d == 1) {
    if (scalar_estimates.empty()) {
      result(0, 0) = std::numeric_limits<double>::quiet_NaN();
      return result;
    }
    std::sort(scalar_estimates.begin(), scalar_estimates.end());
    auto const quantile_type7 = [&scalar_estimates](double probability) {
      double const index =
          probability * static_cast<double>(scalar_estimates.size() - 1);
      std::size_t const lower = static_cast<std::size_t>(std::floor(index));
      std::size_t const upper = static_cast<std::size_t>(std::ceil(index));
      if (lower == upper) return scalar_estimates[lower];
      double const fraction = index - static_cast<double>(lower);
      return scalar_estimates[lower] +
             fraction * (scalar_estimates[upper] - scalar_estimates[lower]);
    };
    double const q25 = quantile_type7(0.25);
    double const q75 = quantile_type7(0.75);
    // variance.lm() defaults outlier_iqr to Inf and still evaluates the
    // threshold. Preserve R's Inf * 0 = NaN behavior for zero-IQR estimates.
    double const threshold = q75 + outlier_iqr * (q75 - q25);
    double total = 0.0;
    std::size_t retained_count = 0;
    for (double const value : scalar_estimates) {
      if (value < threshold) {
        total += value;
        ++retained_count;
      }
    }
    result(0, 0) = retained_count == 0
                       ? std::numeric_limits<double>::quiet_NaN()
                       : total / static_cast<double>(retained_count);
    return result;
  }
  for (unsigned int j = 0; j < d; ++j) {
    for (unsigned int k = 0; k < d; ++k) {
      result(j, k) = counts(j, k) == 0
                         ? std::numeric_limits<double>::quiet_NaN()
                         : sums(j, k) / static_cast<double>(counts(j, k));
    }
  }
  return result;
}

arma::mat estimate_mgaussian_variance(arma::mat const& data,
                                      unsigned int p_response) {
  unsigned int const q = p_response > 0 ? p_response : data.n_cols;
  if (q == 0 || q > data.n_cols) {
    throw std::invalid_argument(
        "fastcpd: p_response must be between 1 and data.n_cols");
  }
  return estimate_linear_regression_variance(data, q);
}

arma::mat resolve_variance_estimate(Options const& options,
                                    arma::mat const& data,
                                    std::string const& family) {
  if (!options.variance_estimate.is_empty()) {
    // Preserve a user-provided covariance exactly when it is well
    // conditioned.  R does not call nearest_pd_() for explicit
    // variance_estimation; it only applies the rcond fallback below.
    arma::mat value = options.variance_estimate;
    if (family == "mean" && value.n_rows != data.n_cols) {
      throw std::invalid_argument(
          "fastcpd: mean variance_estimate must be data.n_cols by data.n_cols");
    }
    if (family == "mgaussian") {
      unsigned int const q =
          options.p_response > 0 ? options.p_response : data.n_cols;
      if (value.n_rows != q) {
        throw std::invalid_argument(
            "fastcpd: mgaussian variance_estimate must be p_response by "
            "p_response");
      }
    }
    return apply_rcond_fallback(std::move(value));
  }
  if (family == "mean") {
    return apply_rcond_fallback(estimate_mean_variance(data));
  }
  if (family == "gaussian") {
    // Gaussian costs do not consume variance_estimate. R only uses this
    // automatic estimate to scale character BIC/MBIC/MDL penalties, so a
    // numeric beta can skip the block-lagged estimator entirely.
    if (options.beta.has_value()) return arma::eye<arma::mat>(1, 1);
    // R's lm/ar wrappers estimate the scalar error variance with
    // variance.lm(data_).  This value is also used to scale a character
    // BIC/MBIC/MDL beta before dispatching the Gaussian SEN solver.
    arma::mat estimated = estimate_linear_regression_variance(data, 1);
    return apply_rcond_fallback(std::move(estimated));
  }
  if (family == "mgaussian") {
    // R's VAR/multivariate-lm path uses variance.lm(), projects the result to
    // the nearest positive-definite matrix, then applies the common rcond
    // guard.  Keep that ordering here as well.
    arma::mat estimated =
        estimate_mgaussian_variance(data, options.p_response);
    if (estimated.is_empty() || !estimated.is_finite()) {
      throw std::invalid_argument(
          "fastcpd: automatic variance estimate is not finite");
    }
    estimated = nearest_positive_definite(estimated);
    return apply_rcond_fallback(std::move(estimated));
  }
  return arma::eye(1, 1);
}

void validate_custom_options(Options const& options,
                             std::string const& family) {
  if (family != "custom") return;
  if (!options.cost_pelt && !options.cost_sen) {
    throw std::invalid_argument(
        "fastcpd: custom family requires cost_pelt or cost_sen");
  }
  if (static_cast<bool>(options.cost_gradient) !=
      static_cast<bool>(options.cost_hessian)) {
    throw std::invalid_argument(
        "fastcpd: provide both cost_gradient and cost_hessian, or neither");
  }
  if (!options.cost_pelt &&
      !(options.cost_sen && options.cost_gradient && options.cost_hessian)) {
    throw std::invalid_argument(
        "fastcpd: custom SEN costs require cost_sen, cost_gradient, and "
        "cost_hessian when cost_pelt is absent");
  }
}

#define FASTCPD_FWD_ARGS                                                     \
  beta, options.cost_pelt, options.cost_sen, options.cost_gradient,           \
      options.cost_hessian, options.cp_only, data, options.epsilon, family,   \
      options.multiple_epochs, line_search, lower, options.momentum_coef,     \
      order, p, p_response, pruning_coef, segment_count, options.trim, upper, \
      vanilla_percentage, variance_estimate, options.warm_start

#define FASTCPD_MAKE_AND_RUN(kRProgress, Policy, kVanillaOnly, kCostAdj,      \
                             kLineSearch, kNDimsValue)                       \
  {                                                                           \
    fastcpd::classes::Fastcpd<                                                \
        fastcpd::families::Policy, kRProgress, kVanillaOnly,                  \
        fastcpd::classes::CostAdjustment::kCostAdj, kLineSearch,              \
        kNDimsValue>                                                          \
        solver(FASTCPD_FWD_ARGS);                                             \
    return solver.Run();                                                      \
  }

#define FASTCPD_DISPATCH_PELT_COST(kRProgress, Policy, kNDimsValue)           \
  if (cost_adjustment == "MBIC") {                                            \
    FASTCPD_MAKE_AND_RUN(kRProgress, Policy, true, kMBIC, false,              \
                         kNDimsValue);                                        \
  }                                                                           \
  if (cost_adjustment == "MDL") {                                             \
    FASTCPD_MAKE_AND_RUN(kRProgress, Policy, true, kMDL, false,               \
                         kNDimsValue);                                        \
  }                                                                           \
  FASTCPD_MAKE_AND_RUN(kRProgress, Policy, true, kBIC, false, kNDimsValue);

#define FASTCPD_DISPATCH_PELT(kRProgress, Policy)                             \
  if (use_1d) {                                                               \
    FASTCPD_DISPATCH_PELT_COST(kRProgress, Policy, 1);                        \
  } else {                                                                    \
    FASTCPD_DISPATCH_PELT_COST(kRProgress, Policy, -1);                       \
  }

#define FASTCPD_DISPATCH_LINE_SEARCH(kRProgress, Policy, kVanillaOnly,        \
                                     kCostAdj)                                \
  if (use_line_search) {                                                       \
    FASTCPD_MAKE_AND_RUN(kRProgress, Policy, kVanillaOnly, kCostAdj, true,    \
                         -1);                                                 \
  }                                                                           \
  FASTCPD_MAKE_AND_RUN(kRProgress, Policy, kVanillaOnly, kCostAdj, false, -1);

#define FASTCPD_DISPATCH_VANILLA(kRProgress, Policy, kCostAdj)                \
  if (vanilla_only) {                                                          \
    FASTCPD_DISPATCH_LINE_SEARCH(kRProgress, Policy, true, kCostAdj);         \
  }                                                                           \
  FASTCPD_DISPATCH_LINE_SEARCH(kRProgress, Policy, false, kCostAdj);

#define FASTCPD_DISPATCH_SEGD(kRProgress, Policy)                             \
  if (cost_adjustment == "MBIC") {                                            \
    FASTCPD_DISPATCH_VANILLA(kRProgress, Policy, kMBIC);                      \
  }                                                                           \
  if (cost_adjustment == "MDL") {                                             \
    FASTCPD_DISPATCH_VANILLA(kRProgress, Policy, kMDL);                       \
  }                                                                           \
  FASTCPD_DISPATCH_VANILLA(kRProgress, Policy, kBIC);

RunResult dispatch(double beta, std::string const& cost_adjustment,
                   Options const& options, arma::mat const& data,
                   std::string const& family, arma::colvec const& order,
                   int p, unsigned int p_response, double pruning_coef,
                   int segment_count, double vanilla_percentage,
                   arma::mat const& variance_estimate,
                   arma::colvec const& lower, arma::colvec const& upper,
                   arma::colvec const& line_search) {
  bool const use_line_search =
      line_search.n_elem > 1 ||
      (line_search.n_elem == 1 && line_search(0) != 1.0);
  bool const vanilla_only = vanilla_percentage == 1.0;
  bool const use_1d = data.n_cols == 1;

  if (family == "mean") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_PELT(true, MeanFamily);
    }
    FASTCPD_DISPATCH_PELT(false, MeanFamily);
  }
  if (family == "variance") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_PELT(true, VarianceFamily);
    }
    FASTCPD_DISPATCH_PELT(false, VarianceFamily);
  }
  if (family == "meanvariance") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_PELT(true, MeanvarianceFamily);
    }
    FASTCPD_DISPATCH_PELT(false, MeanvarianceFamily);
  }
  if (family == "exponential") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_PELT_COST(true, ExponentialFamily, -1);
    }
    FASTCPD_DISPATCH_PELT_COST(false, ExponentialFamily, -1);
  }
  if (family == "mgaussian") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_PELT_COST(true, MgaussianFamily, -1);
    }
    FASTCPD_DISPATCH_PELT_COST(false, MgaussianFamily, -1);
  }
  if (family == "garch") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_PELT_COST(true, GarchFamily, -1);
    }
    FASTCPD_DISPATCH_PELT_COST(false, GarchFamily, -1);
  }
  if (family == "arima") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_PELT_COST(true, ArimaFamily, -1);
    }
    FASTCPD_DISPATCH_PELT_COST(false, ArimaFamily, -1);
  }
  if (family == "gaussian") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_SEGD(true, GaussianFamily);
    }
    FASTCPD_DISPATCH_SEGD(false, GaussianFamily);
  }
  if (family == "lasso") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_SEGD(true, LassoFamily);
    }
    FASTCPD_DISPATCH_SEGD(false, LassoFamily);
  }
  if (family == "binomial") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_SEGD(true, BinomialFamily);
    }
    FASTCPD_DISPATCH_SEGD(false, BinomialFamily);
  }
  if (family == "poisson") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_SEGD(true, PoissonFamily);
    }
    FASTCPD_DISPATCH_SEGD(false, PoissonFamily);
  }
  if (family == "arma") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_SEGD(true, ArmaFamily);
    }
    FASTCPD_DISPATCH_SEGD(false, ArmaFamily);
  }
  if (family == "ma") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_SEGD(true, MaFamily);
    }
    FASTCPD_DISPATCH_SEGD(false, MaFamily);
  }
  if (family == "quantile") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_SEGD(true, QuantileFamily);
    }
    FASTCPD_DISPATCH_SEGD(false, QuantileFamily);
  }
  if (family == "custom") {
    if (options.show_progress) {
      FASTCPD_DISPATCH_SEGD(true, CustomFamily);
    }
    FASTCPD_DISPATCH_SEGD(false, CustomFamily);
  }
  throw std::invalid_argument("fastcpd: unsupported family '" + family + "'");
}

#undef FASTCPD_FWD_ARGS
#undef FASTCPD_MAKE_AND_RUN
#undef FASTCPD_DISPATCH_PELT_COST
#undef FASTCPD_DISPATCH_PELT
#undef FASTCPD_DISPATCH_LINE_SEARCH
#undef FASTCPD_DISPATCH_VANILLA
#undef FASTCPD_DISPATCH_SEGD

Result make_result(RunResult&& value) {
  return Result{std::move(std::get<0>(value)), std::move(std::get<1>(value)),
                std::move(std::get<2>(value)), std::move(std::get<3>(value)),
                std::move(std::get<4>(value))};
}

}  // namespace detail

namespace {

class RRandom {
 public:
  explicit RRandom(std::int32_t seed) {
    if (seed == std::numeric_limits<std::int32_t>::min()) {
      throw std::invalid_argument(
          "fastcpd: R-compatible seed must be greater than -2147483648");
    }
    std::uint32_t value = static_cast<std::uint32_t>(seed);
    for (int index = 0; index < 50; ++index) value = 69069u * value + 1u;
    for (std::size_t index = 0; index <= state_.size(); ++index) {
      value = 69069u * value + 1u;
      if (index > 0) state_[index - 1] = value;
    }
  }

  double uniform() {
    double value = static_cast<double>(raw()) * 0x1p-32;
    if (value == 0.0) {
      value = 0.5 / static_cast<double>(
                        std::numeric_limits<std::uint32_t>::max());
    }
    return value;
  }

  double normal() {
    constexpr double precision = 134217728.0;
    double const probability =
        (std::floor(precision * uniform()) + uniform()) / precision;
    return standard_normal_quantile(probability);
  }

  std::vector<arma::uword> sample_without_replacement(arma::uword population,
                                                       arma::uword size) {
    if (size > population) {
      throw std::invalid_argument(
          "fastcpd: sample size exceeds the KCP bandwidth population");
    }
    std::vector<arma::uword> available(population);
    std::iota(available.begin(), available.end(), 0u);
    std::vector<arma::uword> result(size);
    arma::uword remaining = population;
    for (arma::uword index = 0; index < size; ++index) {
      arma::uword const selected = uniform_index(remaining);
      result[index] = available[selected];
      --remaining;
      available[selected] = available[remaining];
    }
    return result;
  }

  static double standard_normal_quantile(double probability) {
    double const q = probability - 0.5;
    if (std::abs(q) <= 0.425) {
      double const r = 0.180625 - q * q;
      double const numerator =
          (((((((2509.0809287301226727 * r + 33430.575583588128105) * r +
                 67265.770927008700853) *
                    r +
                45921.953931549871457) *
                   r +
               13731.693765509461125) *
                  r +
              1971.5909503065514427) *
                 r +
             133.14166789178437745) *
                r +
            3.387132872796366608);
      double const denominator =
          (((((((5226.495278852854561 * r + 28729.085735721942674) * r +
                 39307.89580009271061) *
                    r +
                21213.794301586595867) *
                   r +
               5394.1960214247511077) *
                  r +
              687.1870074920579083) *
                 r +
             42.313330701600911252) *
                r +
            1.0);
      return q * numerator / denominator;
    }

    double r = q > 0.0 ? 1.0 - probability : probability;
    r = std::sqrt(-std::log(r));
    double value;
    if (r <= 5.0) {
      r -= 1.6;
      double const numerator =
          (((((((7.7454501427834140764e-4 * r +
                 0.0227238449892691845833) *
                    r +
                0.24178072517745061177) *
                   r +
               1.27045825245236838258) *
                  r +
              3.64784832476320460504) *
                 r +
             5.7694972214606914055) *
                r +
            4.6303378461565452959) *
               r +
           1.42343711074968357734);
      double const denominator =
          (((((((1.05075007164441684324e-9 * r +
                 5.475938084995344946e-4) *
                    r +
                0.0151986665636164571966) *
                   r +
               0.14810397642748007459) *
                  r +
              0.68976733498510000455) *
                 r +
             1.6763848301838038494) *
                r +
            2.05319162663775882187) *
               r +
           1.0);
      value = numerator / denominator;
    } else {
      r -= 5.0;
      double const numerator =
          (((((((2.01033439929228813265e-7 * r +
                 2.71155556874348757815e-5) *
                    r +
                0.0012426609473880784386) *
                   r +
               0.026532189526576123093) *
                  r +
              0.29656057182850489123) *
                 r +
             1.7848265399172913358) *
                r +
            5.4637849111641143699) *
               r +
           6.6579046435011037772);
      double const denominator =
          (((((((2.04426310338993978564e-15 * r +
                 1.4215117583164458887e-7) *
                    r +
                1.8463183175100546818e-5) *
                   r +
               7.868691311456132591e-4) *
                  r +
              0.0148753612908506148525) *
                 r +
             0.13692988092273580531) *
                r +
            0.59983220655588793769) *
               r +
           1.0);
      value = numerator / denominator;
    }
    return q < 0.0 ? -value : value;
  }

  std::vector<arma::uword> sample_with_replacement(arma::uword population,
                                                    arma::uword size) {
    if (population == 0 && size > 0) {
      throw std::invalid_argument("fastcpd: cannot sample an empty segment");
    }
    std::vector<arma::uword> result(size);
    for (arma::uword index = 0; index < size; ++index) {
      result[index] = uniform_index(population);
    }
    return result;
  }

 private:
  std::uint32_t raw() {
    if (position_ == state_.size()) twist();
    std::uint32_t value = state_[position_++];
    value ^= value >> 11;
    value ^= (value << 7) & 0x9d2c5680u;
    value ^= (value << 15) & 0xefc60000u;
    value ^= value >> 18;
    return value;
  }

  void twist() {
    constexpr std::uint32_t upper_mask = 0x80000000u;
    constexpr std::uint32_t lower_mask = 0x7fffffffu;
    for (std::size_t index = 0; index < state_.size(); ++index) {
      std::uint32_t const combined =
          (state_[index] & upper_mask) |
          (state_[(index + 1) % state_.size()] & lower_mask);
      state_[index] = state_[(index + 397) % state_.size()] ^ (combined >> 1);
      if (combined & 1u) state_[index] ^= 0x9908b0dfu;
    }
    position_ = 0;
  }

  arma::uword uniform_index(arma::uword population) {
    unsigned int bits = 0;
    for (arma::uword value = population - 1; value > 0; value >>= 1) ++bits;
    unsigned int const chunks = bits / 16 + 1;
    std::uint64_t const mask =
        bits == 0 ? 0u : (std::uint64_t{1} << bits) - 1u;
    while (true) {
      std::uint64_t candidate = 0;
      for (unsigned int chunk = 0; chunk < chunks; ++chunk) {
        candidate = (candidate << 16) + (raw() >> 16);
      }
      candidate &= mask;
      if (candidate < population) return static_cast<arma::uword>(candidate);
    }
  }

  std::array<std::uint32_t, 624> state_{};
  std::size_t position_ = 624;
};

struct KernelSpec {
  arma::uword feature_count;
  double bandwidth;
  arma::colvec public_order;
};

KernelSpec kernel_spec(arma::colvec const& order) {
  double const feature_value = order.n_elem >= 1 ? order(0) : 0.0;
  double const bandwidth = order.n_elem >= 2 ? order(1) : 0.0;
  if (!std::isfinite(feature_value) ||
      (feature_value > 0.0 && feature_value != std::floor(feature_value))) {
    throw std::invalid_argument(
        "fastcpd: KCP feature count must be a positive integer when positive");
  }
  if (!std::isfinite(bandwidth)) {
    throw std::invalid_argument("fastcpd: KCP bandwidth must be finite");
  }
  double const normalized_features = feature_value > 0.0 ? feature_value : 100.0;
  if (normalized_features >
      static_cast<double>(std::numeric_limits<std::int32_t>::max())) {
    throw std::invalid_argument("fastcpd: KCP feature count is too large");
  }
  return KernelSpec{static_cast<arma::uword>(normalized_features), bandwidth,
                    arma::colvec{normalized_features, bandwidth}};
}

double median(std::vector<double>* values) {
  std::sort(values->begin(), values->end());
  std::size_t const middle = values->size() / 2;
  if (values->size() % 2 == 1) return (*values)[middle];
  return ((*values)[middle - 1] + (*values)[middle]) / 2.0;
}

double kernel_bandwidth(arma::mat const& data, double supplied,
                        RRandom* random) {
  if (supplied > 0.0) return supplied;
  arma::mat sampled;
  if (data.n_rows > 1000) {
    std::vector<arma::uword> const indices =
        random->sample_without_replacement(data.n_rows, 1000);
    arma::uvec selected(indices.size());
    for (arma::uword index = 0; index < selected.n_elem; ++index) {
      selected(index) = indices[index];
    }
    sampled = data.rows(selected);
  } else {
    sampled = data;
  }

  std::vector<double> positive_distances;
  positive_distances.reserve(
      static_cast<std::size_t>(sampled.n_rows * (sampled.n_rows - 1) / 2));
  for (arma::uword row = 0; row < sampled.n_rows; ++row) {
    for (arma::uword other = row + 1; other < sampled.n_rows; ++other) {
      double const distance =
          arma::accu(arma::square(sampled.row(row) - sampled.row(other)));
      if (distance > 0.0) positive_distances.push_back(distance);
    }
  }
  return positive_distances.empty()
             ? 1.0
             : std::sqrt(median(&positive_distances) / 2.0);
}

arma::mat kernel_transform(arma::mat const& data, KernelSpec const& spec,
                           RRandom* random) {
  double const bandwidth = kernel_bandwidth(data, spec.bandwidth, random);
  arma::mat omega(data.n_cols, spec.feature_count);
  for (arma::uword index = 0; index < omega.n_elem; ++index) {
    omega(index) = random->normal() / bandwidth;
  }
  arma::rowvec phase(spec.feature_count);
  for (arma::uword index = 0; index < phase.n_elem; ++index) {
    phase(index) = 2.0 * arma::datum::pi * random->uniform();
  }
  arma::mat projection = data * omega;
  projection.each_row() += phase;
  return std::sqrt(2.0 / static_cast<double>(spec.feature_count)) *
         arma::cos(projection);
}

std::int32_t kernel_seed(std::optional<std::int32_t> seed) {
  if (seed.has_value()) return *seed;
  std::random_device random_device;
  std::uint32_t const value = random_device() & 0x7fffffffu;
  return static_cast<std::int32_t>(value);
}

Result detect_native(arma::mat const& data, Options options) {
  if (data.n_rows == 0 || data.n_cols == 0) {
    throw std::invalid_argument("fastcpd: data must be a non-empty matrix");
  }
  if (options.cost_adjustment.empty()) options.cost_adjustment = "MBIC";
  detail::validate_common_options(options);
  detail::validate_cost_adjustment(options.cost_adjustment);
  std::string const family =
      detail::normalize_family(options.family, options.order);
  detail::validate_custom_options(options, family);

  arma::colvec const order = options.order;
  detail::validate_family_input(family, data, order, options.p_response);
  int const p = detail::compute_p(family, data, order, options.p_response,
                                  options.p);
  if (p <= 0) {
    throw std::invalid_argument(
        "fastcpd: inferred p must be positive; pass Options::p explicitly");
  }

  unsigned int const p_response =
      family == "mgaussian" && options.p_response == 0
          ? static_cast<unsigned int>(data.n_cols)
          : options.p_response;
  arma::colvec const line_search =
      detail::normalize_line_search(options.line_search);
  arma::colvec const lower = detail::fill_or_validate_bound(
      options.lower, p, -std::numeric_limits<double>::infinity(), "lower");
  arma::colvec const upper = detail::fill_or_validate_bound(
      options.upper, p, std::numeric_limits<double>::infinity(), "upper");
  if (arma::any(lower > upper)) {
    throw std::invalid_argument(
        "fastcpd: lower values must not exceed upper values");
  }
  arma::mat const variance_estimate =
      detail::resolve_variance_estimate(options, data, family);

  double beta = detail::compute_beta(
      options.beta, options.beta_criterion, static_cast<int>(data.n_rows), p);
  // R scales character BIC/MBIC/MDL penalties by the estimated Gaussian
  // variance.  Numeric beta values are already on the caller's scale and are
  // intentionally left untouched.  `family == "gaussian"` also covers the
  // R/Python AR and lm wrappers after their design transformation.
  if (!options.beta.has_value() && family == "gaussian") {
    if (variance_estimate.n_elem != 1) {
      throw std::invalid_argument(
          "fastcpd: variance_estimate for Gaussian criterion beta must be "
          "scalar");
    }
    beta *= variance_estimate(0, 0);
  }
  double const pruning_coef = detail::compute_pruning_coef(
      options.pruning_coef, options.cost_adjustment, family, p);
  double const vanilla_percentage =
      detail::is_pelt_family(family) ? 1.0 : options.vanilla_percentage;

  Result result = detail::make_result(detail::dispatch(
      beta, options.cost_adjustment, options, data, family, order, p,
      p_response, pruning_coef, options.segment_count, vanilla_percentage,
      variance_estimate, lower, upper, line_search));
  result.family = family;
  result.order = order;
  result.cp_only = options.cp_only;
  return result;
}

Result with_public_metadata(Result result, std::string family,
                            arma::colvec order, bool cp_only) {
  result.family = std::move(family);
  result.order = std::move(order);
  result.cp_only = cp_only;
  return result;
}

void require_finite_data(arma::mat const& data) {
  if (data.n_rows == 0 || data.n_cols == 0) {
    throw std::invalid_argument("fastcpd: data must be a non-empty matrix");
  }
  if (!data.is_finite()) {
    throw std::invalid_argument(
        "fastcpd: data must contain only finite values");
  }
}

void require_univariate(arma::mat const& data, char const* family) {
  if (data.n_cols != 1) {
    throw std::invalid_argument(std::string("fastcpd: ") + family +
                                " data must be univariate");
  }
}

unsigned int positive_scalar_order(arma::colvec const& order,
                                   char const* family) {
  if (order.n_elem != 1 || !detail::is_nonnegative_integer(order(0)) ||
      order(0) <= 0.0) {
    throw std::invalid_argument(std::string("fastcpd: ") + family +
                                " order must be a positive integer");
  }
  return static_cast<unsigned int>(order(0));
}

unsigned int ar_order(arma::colvec const& order) {
  if (order.n_elem == 1) return positive_scalar_order(order, "AR");
  detail::validate_integer_order(order, 3, "AR");
  if (order(0) <= 0.0 || order(1) != 0.0 || order(2) != 0.0) {
    throw std::invalid_argument(
        "fastcpd: AR order must be p or (p, 0, 0), with p positive");
  }
  return static_cast<unsigned int>(order(0));
}

arma::mat make_ar_design(arma::colvec const& data, unsigned int order) {
  if (data.n_rows <= order) {
    throw std::invalid_argument(
        "fastcpd: AR order must be smaller than data.n_rows");
  }
  arma::uword const rows = data.n_rows - order;
  arma::mat design(rows, order + 1);
  design.col(0) = data.rows(order, data.n_rows - 1);
  for (unsigned int lag = 1; lag <= order; ++lag) {
    design.col(lag) = data.rows(order - lag, data.n_rows - lag - 1);
  }
  return design;
}

arma::mat make_var_design(arma::mat const& data, unsigned int order) {
  if (data.n_rows <= order) {
    throw std::invalid_argument(
        "fastcpd: VAR order must be smaller than data.n_rows");
  }
  arma::uword const rows = data.n_rows - order;
  arma::uword const q = data.n_cols;
  arma::mat design(rows, q * (order + 1));
  design.cols(0, q - 1) = data.rows(order, data.n_rows - 1);
  for (unsigned int lag = 1; lag <= order; ++lag) {
    arma::uword const first = q * lag;
    design.cols(first, first + q - 1) =
        data.rows(order - lag, data.n_rows - lag - 1);
  }
  return design;
}

Result restore_lag_coordinates(Result result, unsigned int offset,
                               arma::uword original_rows,
                               arma::uword response_count) {
  result.raw_change_points += static_cast<double>(offset);
  result.change_points += static_cast<double>(offset);
  if (!result.residuals.is_empty()) {
    arma::mat residuals(
        original_rows, response_count,
        arma::fill::value(std::numeric_limits<double>::quiet_NaN()));
    residuals.rows(offset, original_rows - 1) = result.residuals;
    result.residuals = std::move(residuals);
  }
  return result;
}

arma::mat centered_ranks(arma::mat const& data) {
  arma::mat result(data.n_rows, data.n_cols);
  std::vector<arma::uword> indices(data.n_rows);
  for (arma::uword column = 0; column < data.n_cols; ++column) {
    for (arma::uword row = 0; row < data.n_rows; ++row) indices[row] = row;
    std::stable_sort(indices.begin(), indices.end(),
                     [&data, column](arma::uword lhs, arma::uword rhs) {
                       return data(lhs, column) < data(rhs, column);
                     });
    arma::uword start = 0;
    while (start < indices.size()) {
      arma::uword end = start + 1;
      while (end < indices.size() &&
             data(indices[end], column) == data(indices[start], column)) {
        ++end;
      }
      double const average_rank =
          (static_cast<double>(start + 1) + static_cast<double>(end)) / 2.0;
      for (arma::uword position = start; position < end; ++position) {
        result(indices[position], column) = average_rank;
      }
      start = end;
    }
  }
  result -= static_cast<double>(data.n_rows + 1) / 2.0;
  return result;
}

Result detect_kernel_with_random(arma::mat const& data, Options options,
                                 RRandom* random) {
  require_finite_data(data);
  KernelSpec const spec = kernel_spec(options.order);
  bool const cp_only = options.cp_only;
  arma::mat transformed = kernel_transform(data, spec, random);
  if (!options.beta.has_value()) {
    options.beta =
        (static_cast<double>(data.n_cols) + 2.0) *
        std::log(static_cast<double>(data.n_rows)) / 2.0;
  }
  if (options.cost_adjustment.empty()) options.cost_adjustment = "BIC";
  options.family = "mean";
  options.order = arma::colvec{0.0, 0.0, 0.0};
  options.p = static_cast<int>(spec.feature_count);
  options.vanilla_percentage = 1.0;
  if (options.variance_estimate.is_empty()) {
    options.variance_estimate =
        arma::eye<arma::mat>(spec.feature_count, spec.feature_count);
  }
  Result result = detect_native(transformed, std::move(options));
  return with_public_metadata(std::move(result), "kcp", spec.public_order,
                              cp_only);
}

}  // namespace

Result detect(arma::mat const& data, Options options) {
  std::string const family = detail::lower_ascii(options.family);
  if (family == "lm") return detect_lm(data, std::move(options));
  if (family == "var") return detect_var(data, std::move(options));
  if (family == "rank") return detect_rank(data, std::move(options));
  if (family == "kcp") return detect_kcp(data, std::move(options));
  if (family == "ar") {
    require_univariate(data, "AR");
    return detect_ar(data.col(0), std::move(options));
  }
  if (family == "arma") {
    require_univariate(data, "ARMA");
    return detect_arma(data.col(0), std::move(options));
  }
  if (family == "arima") {
    require_univariate(data, "ARIMA");
    return detect_arima(data.col(0), std::move(options));
  }
  if (family == "garch") {
    require_univariate(data, "GARCH");
    return detect_garch(data.col(0), std::move(options));
  }
  return detect_native(data, std::move(options));
}

Result detect(arma::colvec const& data, Options options) {
  return detect(arma::mat(data), std::move(options));
}

Result detect_mean(arma::mat const& data, Options options) {
  bool const cp_only = options.cp_only;
  options.family = "mean";
  options.order = arma::colvec{0.0, 0.0, 0.0};
  return with_public_metadata(detect_native(data, std::move(options)), "mean",
                              arma::colvec{0.0, 0.0, 0.0}, cp_only);
}

Result detect_variance(arma::mat const& data, Options options) {
  bool const cp_only = options.cp_only;
  options.family = "variance";
  options.order = arma::colvec{0.0, 0.0, 0.0};
  return with_public_metadata(detect_native(data, std::move(options)),
                              "variance", arma::colvec{0.0, 0.0, 0.0},
                              cp_only);
}

Result detect_meanvariance(arma::mat const& data, Options options) {
  bool const cp_only = options.cp_only;
  options.family = "meanvariance";
  options.order = arma::colvec{0.0, 0.0, 0.0};
  return with_public_metadata(detect_native(data, std::move(options)),
                              "meanvariance", arma::colvec{0.0, 0.0, 0.0},
                              cp_only);
}

Result detect_mean_variance(arma::mat const& data, Options options) {
  return detect_meanvariance(data, std::move(options));
}

Result detect_exponential(arma::mat const& data, Options options) {
  bool const cp_only = options.cp_only;
  options.family = "exponential";
  options.order = arma::colvec{0.0, 0.0, 0.0};
  return with_public_metadata(detect_native(data, std::move(options)),
                              "exponential", arma::colvec{0.0, 0.0, 0.0},
                              cp_only);
}

Result detect_lm(arma::mat const& data, Options options) {
  if (data.n_cols < 2) {
    throw std::invalid_argument(
        "fastcpd: linear-regression data must include predictors");
  }
  unsigned int const q = options.p_response > 0 ? options.p_response : 1u;
  if (q >= data.n_cols) {
    throw std::invalid_argument(
        "fastcpd: linear-regression response columns must leave predictors");
  }
  bool const cp_only = options.cp_only;
  arma::colvec const public_order = options.order;
  if (q > 1) {
    options.family = "mgaussian";
    options.p_response = q;
    if (options.p == 0) {
      options.p = static_cast<int>(q * (data.n_cols - q));
    }
    options.vanilla_percentage = 1.0;
  } else {
    options.family = "gaussian";
    options.p_response = 1;
    if (options.p == 0) options.p = static_cast<int>(data.n_cols - 1);
  }
  return with_public_metadata(detect_native(data, std::move(options)), "lm",
                              public_order, cp_only);
}

Result detect_linear_regression(arma::mat const& data, Options options) {
  return detect_lm(data, std::move(options));
}

#define FASTCPD_DEFINE_REGRESSION_WRAPPER(Name, Family)                      \
  Result Name(arma::mat const& data, Options options) {                       \
    bool const cp_only = options.cp_only;                                     \
    arma::colvec const public_order = options.order;                          \
    options.family = Family;                                                  \
    return with_public_metadata(detect_native(data, std::move(options)),      \
                                Family, public_order, cp_only);               \
  }

FASTCPD_DEFINE_REGRESSION_WRAPPER(detect_lasso, "lasso")
FASTCPD_DEFINE_REGRESSION_WRAPPER(detect_binomial, "binomial")
FASTCPD_DEFINE_REGRESSION_WRAPPER(detect_poisson, "poisson")
FASTCPD_DEFINE_REGRESSION_WRAPPER(detect_quantile, "quantile")

#undef FASTCPD_DEFINE_REGRESSION_WRAPPER

Result detect_logistic_regression(arma::mat const& data, Options options) {
  return detect_binomial(data, std::move(options));
}

Result detect_poisson_regression(arma::mat const& data, Options options) {
  return detect_poisson(data, std::move(options));
}

Result detect_quantile_regression(arma::mat const& data, Options options) {
  return detect_quantile(data, std::move(options));
}

Result detect_ar(arma::colvec const& data, Options options) {
  arma::colvec const public_order = options.order;
  unsigned int const order = ar_order(public_order);
  bool const cp_only = options.cp_only;
  arma::mat design = make_ar_design(data, order);
  options.family = "gaussian";
  options.order = arma::colvec{static_cast<double>(order)};
  options.p = static_cast<int>(order);
  options.p_response = 1;
  Result result = detect_native(design, std::move(options));
  result = restore_lag_coordinates(std::move(result), order, data.n_rows, 1);
  return with_public_metadata(std::move(result), "ar", public_order, cp_only);
}

Result detect_arma(arma::colvec const& data, Options options) {
  arma::colvec const public_order = options.order;
  detail::validate_integer_order(public_order, 2, "ARMA");
  unsigned int const p = static_cast<unsigned int>(public_order(0));
  unsigned int const q = static_cast<unsigned int>(public_order(1));
  if (p == 0 && q == 0) {
    throw std::invalid_argument(
        "fastcpd: ARMA order must contain a non-zero value");
  }
  bool const cp_only = options.cp_only;
  Result result;
  if (q == 0) {
    options.order = arma::colvec{static_cast<double>(p)};
    result = detect_ar(data, std::move(options));
  } else {
    options.family = p == 0 ? "ma" : "arma";
    options.p = static_cast<int>(p + q + 1);
    result = detect_native(arma::mat(data), std::move(options));
  }
  return with_public_metadata(std::move(result), "arma", public_order,
                              cp_only);
}

Result detect_arima(arma::colvec const& data, Options options) {
  if (options.include_mean) {
    throw std::invalid_argument(
        "fastcpd: include_mean=true is unsupported by the zero-mean ARIMA "
        "likelihood");
  }
  arma::colvec const public_order = options.order;
  detail::validate_integer_order(public_order, 3, "ARIMA");
  if (arma::all(public_order == 0.0)) {
    throw std::invalid_argument(
        "fastcpd: ARIMA order must contain a non-zero value");
  }
  unsigned int const p = static_cast<unsigned int>(public_order(0));
  unsigned int const d = static_cast<unsigned int>(public_order(1));
  unsigned int const q = static_cast<unsigned int>(public_order(2));
  if (d >= data.n_rows) {
    throw std::invalid_argument(
        "fastcpd: ARIMA integration order must be smaller than data.n_rows");
  }
  bool const cp_only = options.cp_only;
  Result result;
  if (d == 0) {
    options.order = arma::colvec{static_cast<double>(p),
                                 static_cast<double>(q)};
    result = detect_arma(data, std::move(options));
  } else {
    options.family = "arima";
    options.p = static_cast<int>(p + q + 1);
    options.vanilla_percentage = 1.0;
    result = detect_native(arma::mat(data), std::move(options));
  }
  return with_public_metadata(std::move(result), "arima", public_order,
                              cp_only);
}

Result detect_garch(arma::colvec const& data, Options options) {
  arma::colvec const public_order = options.order;
  detail::validate_integer_order(public_order, 2, "GARCH");
  if (arma::all(public_order == 0.0)) {
    throw std::invalid_argument(
        "fastcpd: GARCH order must contain a non-zero value");
  }
  bool const cp_only = options.cp_only;
  options.family = "garch";
  options.p = static_cast<int>(arma::sum(public_order)) + 1;
  options.vanilla_percentage = 1.0;
  return with_public_metadata(
      detect_native(arma::mat(data), std::move(options)), "garch",
      public_order, cp_only);
}

Result detect_var(arma::mat const& data, Options options) {
  require_finite_data(data);
  if (options.p_response != 0) {
    throw std::invalid_argument(
        "fastcpd: p_response is not accepted for raw VAR input");
  }
  arma::colvec const public_order = options.order;
  unsigned int const order = positive_scalar_order(public_order, "VAR");
  bool const cp_only = options.cp_only;
  unsigned int const q = static_cast<unsigned int>(data.n_cols);
  arma::mat design = make_var_design(data, order);
  options.family = "mgaussian";
  options.order = arma::colvec{static_cast<double>(order)};
  options.p_response = q;
  options.p = static_cast<int>(order * q * q);
  options.vanilla_percentage = 1.0;
  Result result = detect_native(design, std::move(options));
  result = restore_lag_coordinates(std::move(result), order, data.n_rows, q);
  return with_public_metadata(std::move(result), "var", public_order, cp_only);
}

Result detect_mgaussian(arma::mat const& data, Options options) {
  bool const cp_only = options.cp_only;
  arma::colvec const public_order = options.order;
  options.family = "mgaussian";
  return with_public_metadata(detect_native(data, std::move(options)),
                              "mgaussian", public_order, cp_only);
}

Result detect_rank(arma::mat const& data, Options options) {
  require_finite_data(data);
  bool const cp_only = options.cp_only;
  arma::colvec const public_order = options.order;
  arma::mat transformed = centered_ranks(data);
  options.family = "mean";
  options.order = arma::colvec{0.0, 0.0, 0.0};
  Result result = detect_native(transformed, std::move(options));
  return with_public_metadata(std::move(result), "rank", public_order,
                              cp_only);
}

Result detect_kernel(arma::mat const& data, Options options) {
  RRandom random(kernel_seed(options.seed));
  return detect_kernel_with_random(data, std::move(options), &random);
}

Result detect_kcp(arma::mat const& data, Options options) {
  return detect_kernel(data, std::move(options));
}

Result mean(arma::mat const& data, Options options) {
  return detect_mean(data, std::move(options));
}

Result variance(arma::mat const& data, Options options) {
  return detect_variance(data, std::move(options));
}

Result meanvariance(arma::mat const& data, Options options) {
  return detect_meanvariance(data, std::move(options));
}

Result exponential(arma::mat const& data, Options options) {
  return detect_exponential(data, std::move(options));
}

Result gaussian(arma::mat const& data, Options options) {
  bool const cp_only = options.cp_only;
  arma::colvec const public_order = options.order;
  options.family = "gaussian";
  return with_public_metadata(detect_native(data, std::move(options)),
                              "gaussian", public_order, cp_only);
}

Result lm(arma::mat const& data, Options options) {
  return detect_lm(data, std::move(options));
}
Result lasso(arma::mat const& data, Options options) {
  return detect_lasso(data, std::move(options));
}
Result binomial(arma::mat const& data, Options options) {
  return detect_binomial(data, std::move(options));
}
Result poisson(arma::mat const& data, Options options) {
  return detect_poisson(data, std::move(options));
}
Result quantile(arma::mat const& data, Options options) {
  return detect_quantile(data, std::move(options));
}
Result ar(arma::colvec const& data, Options options) {
  return detect_ar(data, std::move(options));
}
Result arima(arma::colvec const& data, Options options) {
  return detect_arima(data, std::move(options));
}
Result var(arma::mat const& data, Options options) {
  return detect_var(data, std::move(options));
}
Result mgaussian(arma::mat const& data, Options options) {
  return detect_mgaussian(data, std::move(options));
}
Result rank(arma::mat const& data, Options options) {
  return detect_rank(data, std::move(options));
}
Result kernel(arma::mat const& data, Options options) {
  return detect_kernel(data, std::move(options));
}
Result kcp(arma::mat const& data, Options options) {
  return detect_kcp(data, std::move(options));
}

arma::mat estimate_variance_mean(arma::mat const& data) {
  require_finite_data(data);
  if (data.n_rows < 2) {
    return arma::mat(
        data.n_cols, data.n_cols,
        arma::fill::value(std::numeric_limits<double>::quiet_NaN()));
  }
  arma::mat const differences =
      data.rows(1, data.n_rows - 1) - data.rows(0, data.n_rows - 2);
  return differences.t() * differences /
         (2.0 * static_cast<double>(differences.n_rows));
}

double estimate_variance_median(arma::mat const& data) {
  require_finite_data(data);
  if (data.n_elem < 2) return std::numeric_limits<double>::quiet_NaN();
  double absolute_difference_sum = 0.0;
  double const* const values = data.memptr();
  for (arma::uword index = 1; index < data.n_elem; ++index) {
    absolute_difference_sum += std::abs(values[index] - values[index - 1]);
  }
  double const mean_absolute_difference =
      absolute_difference_sum / static_cast<double>(data.n_elem - 1);
  double const scaled = 2.0 * mean_absolute_difference / 3.0;
  return 2.0 * scaled * scaled;
}

arma::mat estimate_variance_linear_regression(arma::mat const& data,
                                               unsigned int d,
                                               unsigned int block_size,
                                               double outlier_iqr) {
  require_finite_data(data);
  if (d == 0 || d > data.n_cols) {
    throw std::invalid_argument(
        "fastcpd: d must be between 1 and data.n_cols");
  }
  if (std::isnan(outlier_iqr)) {
    throw std::invalid_argument("fastcpd: outlier_iqr must not be NaN");
  }
  return detail::estimate_linear_regression_variance(
      data, d, block_size, outlier_iqr);
}

arma::mat estimate_variance_lm(arma::mat const& data, unsigned int d,
                               unsigned int block_size, double outlier_iqr) {
  return estimate_variance_linear_regression(data, d, block_size,
                                              outlier_iqr);
}

VarianceArmaResult estimate_variance_arma(arma::colvec const& data,
                                          unsigned int p, unsigned int q,
                                          unsigned int max_order) {
  if (data.n_elem == 0 || !data.is_finite()) {
    throw std::invalid_argument(
        "fastcpd: data must be non-empty and contain only finite values");
  }
  if (max_order == 0) max_order = p * q;
  if (max_order == 0) {
    throw std::invalid_argument("fastcpd: max_order must be positive");
  }
  if (data.n_elem <= static_cast<arma::uword>(max_order + 1u)) {
    throw std::invalid_argument(
        "fastcpd: data must contain more observations than max_order");
  }

  VarianceArmaResult result;
  result.table.reserve(max_order);
  double best_aic = std::numeric_limits<double>::infinity();
  double best_bic = std::numeric_limits<double>::infinity();
  result.sigma2_aic = std::numeric_limits<double>::quiet_NaN();
  result.sigma2_bic = std::numeric_limits<double>::quiet_NaN();
  for (unsigned int order = 1; order <= max_order; ++order) {
    arma::mat const design = make_ar_design(data, order);
    arma::mat const estimate = estimate_variance_linear_regression(design);
    double const sigma2 = estimate(0, 0);
    double aic = std::numeric_limits<double>::quiet_NaN();
    double bic = std::numeric_limits<double>::quiet_NaN();
    if (std::isfinite(sigma2) && sigma2 > 0.0) {
      aic = std::log(sigma2) +
            2.0 * static_cast<double>(order) /
                static_cast<double>(data.n_elem);
      bic = std::log(sigma2) +
            static_cast<double>(order) *
                std::log(static_cast<double>(data.n_elem)) /
                static_cast<double>(data.n_elem);
      if (aic < best_aic) {
        best_aic = aic;
        result.sigma2_aic = sigma2;
      }
      if (bic < best_bic) {
        best_bic = bic;
        result.sigma2_bic = sigma2;
      }
    }
    result.table.push_back(VarianceArmaRow{order, sigma2, aic, bic});
  }
  return result;
}

namespace {

using ProfileCost = std::function<double(arma::uword, arma::uword)>;

std::vector<arma::uword> confidence_change_points(Result const& result,
                                                  arma::uword data_rows) {
  std::vector<arma::uword> points;
  points.reserve(result.change_points.n_elem);
  for (double const value : result.change_points) {
    if (!std::isfinite(value) || value <= 0.0 ||
        value >= static_cast<double>(data_rows) ||
        value != std::floor(value)) {
      throw std::invalid_argument(
          "fastcpd: result change points must be interior integer boundaries");
    }
    points.push_back(static_cast<arma::uword>(value));
  }
  std::sort(points.begin(), points.end());
  if (std::adjacent_find(points.begin(), points.end()) != points.end()) {
    throw std::invalid_argument(
        "fastcpd: result change points must be unique boundaries");
  }
  return points;
}

arma::mat confidence_analysis_data(Result const& result,
                                   arma::mat const& data) {
  if (result.family == "rank") return centered_ranks(data);
  return data;
}

std::string confidence_native_family(std::string family) {
  family = detail::lower_ascii(std::move(family));
  if (family == "rank") return "mean";
  if (family == "gaussian") return "lm";
  if (family == "mgaussian") return "var";
  return family;
}

double stable_log_determinant(arma::mat const& matrix) {
  if (!matrix.is_finite()) return std::numeric_limits<double>::infinity();
  arma::vec eigenvalues;
  if (!arma::eig_sym(eigenvalues, arma::symmatl(matrix))) {
    return std::numeric_limits<double>::infinity();
  }
  double total = 0.0;
  for (double const value : eigenvalues) {
    total += std::log(
        std::max(value, std::numeric_limits<double>::epsilon()));
  }
  return total;
}

arma::mat confidence_precision(arma::mat const& covariance,
                               arma::uword dimension) {
  if (covariance.n_rows != dimension || covariance.n_cols != dimension ||
      !covariance.is_finite()) {
    return arma::eye<arma::mat>(dimension, dimension);
  }
  arma::mat precision;
  if (!arma::inv(precision, covariance) || !precision.is_finite()) {
    return arma::eye<arma::mat>(dimension, dimension);
  }
  return precision;
}

Result single_segment_fit(arma::mat const& data, std::string const& family,
                          arma::colvec const& order,
                          Options detector_options) {
  detector_options.family = family;
  detector_options.order = order;
  detector_options.beta = 1e100;
  detector_options.cost_adjustment = "BIC";
  detector_options.cp_only = false;
  detector_options.segment_count = 1;
  detector_options.trim = 0.0;
  detector_options.vanilla_percentage = 1.0;
  detector_options.pruning_coef.reset();
  detector_options.show_progress = false;
  return detect(data, std::move(detector_options));
}

double single_segment_cost(arma::mat const& data, std::string const& family,
                           arma::colvec const& order,
                           Options const& detector_options) {
  try {
    Result const fit =
        single_segment_fit(data, family, order, detector_options);
    if (fit.change_points.n_elem != 0 || fit.cost_values.n_elem != 1) {
      return std::numeric_limits<double>::infinity();
    }
    return fit.cost_values(0);
  } catch (...) {
    return std::numeric_limits<double>::infinity();
  }
}

ProfileCost profile_cost_function(Result const& result,
                                  arma::mat const& analysis_data,
                                  Options const& detector_options) {
  std::string const family = confidence_native_family(result.family);
  if (family == "mean") {
    arma::mat covariance = detector_options.variance_estimate;
    if (covariance.is_empty()) covariance = estimate_variance_mean(analysis_data);
    arma::mat const precision =
        confidence_precision(covariance, analysis_data.n_cols);
    return [&analysis_data, precision](arma::uword start, arma::uword end) {
      arma::mat const segment = analysis_data.rows(start, end - 1);
      arma::mat const centered =
          segment.each_row() - arma::mean(segment, 0);
      return arma::accu((centered * precision) % centered) / 2.0;
    };
  }
  if (family == "variance") {
    arma::mat const centered_data =
        analysis_data.each_row() - arma::mean(analysis_data, 0);
    return [centered_data](arma::uword start, arma::uword end) {
      arma::mat const segment = centered_data.rows(start, end - 1);
      arma::mat const covariance = segment.t() * segment /
                                   static_cast<double>(segment.n_rows);
      return static_cast<double>(segment.n_rows) *
             stable_log_determinant(covariance) / 2.0;
    };
  }
  if (family == "meanvariance") {
    return [&analysis_data](arma::uword start, arma::uword end) {
      arma::mat const segment = analysis_data.rows(start, end - 1);
      arma::mat const centered =
          segment.each_row() - arma::mean(segment, 0);
      arma::mat const covariance = centered.t() * centered /
                                   static_cast<double>(segment.n_rows);
      return static_cast<double>(segment.n_rows) *
             stable_log_determinant(covariance) / 2.0;
    };
  }
  if (family == "exponential") {
    return [&analysis_data](arma::uword start, arma::uword end) {
      arma::colvec const segment = analysis_data.col(0).rows(start, end - 1);
      if (arma::any(segment <= 0.0)) {
        return std::numeric_limits<double>::infinity();
      }
      return static_cast<double>(segment.n_elem) *
             (std::log(arma::mean(segment)) + 1.0);
    };
  }
  if (family == "lm") {
    double sigma2 = 1.0;
    if (!detector_options.variance_estimate.is_empty()) {
      sigma2 = detector_options.variance_estimate(0, 0);
    } else {
      try {
        sigma2 = estimate_variance_linear_regression(analysis_data)(0, 0);
      } catch (...) {
        sigma2 = 1.0;
      }
    }
    if (!std::isfinite(sigma2) || sigma2 <= 0.0) sigma2 = 1.0;
    return [&analysis_data, sigma2](arma::uword start, arma::uword end) {
      arma::mat const segment = analysis_data.rows(start, end - 1);
      arma::colvec const y = segment.col(0);
      arma::mat const x = segment.cols(1, segment.n_cols - 1);
      if (x.n_rows <= x.n_cols) {
        return std::numeric_limits<double>::infinity();
      }
      arma::colvec coefficients;
      if (!arma::solve(coefficients, x, y) || !coefficients.is_finite()) {
        return std::numeric_limits<double>::infinity();
      }
      arma::colvec const residuals = y - x * coefficients;
      return arma::dot(residuals, residuals) / (2.0 * sigma2);
    };
  }
  if (family == "binomial" || family == "poisson" ||
      family == "quantile" || family == "arima") {
    arma::colvec const order = result.order;
    return [&analysis_data, family, order, detector_options](
               arma::uword start, arma::uword end) {
      return single_segment_cost(analysis_data.rows(start, end - 1), family,
                                 order, detector_options);
    };
  }
  throw std::invalid_argument("fastcpd: profile intervals are unavailable for " +
                              result.family);
}

std::vector<ConfidenceInterval> confidence_profile(
    Result const& result, arma::mat const& data,
    ConfidenceOptions const& options) {
  std::vector<arma::uword> const points =
      confidence_change_points(result, data.n_rows);
  if (points.empty()) return {};
  arma::mat const analysis_data = confidence_analysis_data(result, data);
  ProfileCost const cost =
      profile_cost_function(result, analysis_data, options.detector_options);
  double const z =
      RRandom::standard_normal_quantile((1.0 + options.level) / 2.0);
  double const cutoff = z * z / 2.0;
  std::vector<arma::uword> bounds;
  bounds.reserve(points.size() + 2);
  bounds.push_back(0);
  bounds.insert(bounds.end(), points.begin(), points.end());
  bounds.push_back(data.n_rows);

  std::vector<ConfidenceInterval> rows;
  rows.reserve(points.size());
  for (std::size_t index = 0; index < points.size(); ++index) {
    arma::uword const point = points[index];
    arma::uword const left = bounds[index];
    arma::uword const right = bounds[index + 2];
    bool const has_candidate_range =
        options.min_segment_length <= (right - left) / 2;
    arma::uword minimum = 0;
    arma::uword maximum = 0;
    if (has_candidate_range) {
      minimum = left + options.min_segment_length;
      maximum = right - options.min_segment_length;
    }
    if (has_candidate_range && options.window.has_value()) {
      arma::uword const window = *options.window;
      minimum = std::max(minimum, point > window ? point - window : 0u);
      maximum = std::min(maximum, point + window);
    }

    ConfidenceInterval row;
    row.parm = "cp";
    row.method = "profile";
    row.index = static_cast<unsigned int>(index + 1);
    row.estimate = static_cast<double>(point);
    row.cutoff = cutoff;
    row.level = options.level;
    if (has_candidate_range && minimum <= maximum) {
      double profile_min = std::numeric_limits<double>::infinity();
      std::vector<std::pair<arma::uword, double>> finite;
      for (arma::uword candidate = minimum; candidate <= maximum; ++candidate) {
        double const value =
            cost(left, candidate) + cost(candidate, right);
        if (std::isfinite(value)) {
          finite.emplace_back(candidate, value);
          profile_min = std::min(profile_min, value);
        }
      }
      if (!finite.empty()) {
        arma::uword lower = point;
        arma::uword upper = point;
        for (auto const& candidate : finite) {
          if (candidate.second - profile_min <= cutoff) {
            lower = std::min(lower, candidate.first);
            upper = std::max(upper, candidate.first);
          }
        }
        row.lower = static_cast<double>(lower);
        row.upper = static_cast<double>(upper);
        row.profile_min = profile_min;
      }
    }
    rows.push_back(std::move(row));
  }
  return rows;
}

arma::colvec missing_standard_errors(arma::uword count) {
  return arma::colvec(
      count, arma::fill::value(std::numeric_limits<double>::quiet_NaN()));
}

arma::colvec theta_standard_errors(arma::mat const& data,
                                   std::string const& family,
                                   arma::colvec const& order,
                                   Options const& detector_options) {
  if (family == "mean") {
    if (data.n_rows <= 1) return missing_standard_errors(data.n_cols);
    arma::mat const centered = data.each_row() - arma::mean(data, 0);
    arma::mat const covariance =
        centered.t() * centered / static_cast<double>(data.n_rows - 1);
    return arma::sqrt(covariance.diag() /
                      static_cast<double>(data.n_rows));
  }
  if (family == "exponential") {
    double const rate = 1.0 / arma::mean(data.col(0));
    return arma::colvec{rate / std::sqrt(static_cast<double>(data.n_rows))};
  }
  if (family == "lm") {
    arma::colvec const y = data.col(0);
    arma::mat const x = data.cols(1, data.n_cols - 1);
    arma::colvec coefficients;
    arma::mat inverse;
    if (!arma::solve(coefficients, x, y) ||
        !arma::inv(inverse, x.t() * x)) {
      return missing_standard_errors(x.n_cols);
    }
    arma::colvec const residuals = y - x * coefficients;
    arma::uword const degrees =
        std::max<arma::uword>(x.n_rows - x.n_cols, 1);
    double const sigma2 =
        arma::dot(residuals, residuals) / static_cast<double>(degrees);
    return arma::sqrt(inverse.diag() * sigma2);
  }
  if (family == "binomial" || family == "poisson") {
    arma::mat const x = data.cols(1, data.n_cols - 1);
    Result fitted;
    try {
      fitted = single_segment_fit(data, family, order, detector_options);
    } catch (...) {
      return missing_standard_errors(x.n_cols);
    }
    if (fitted.thetas.n_cols != 1 || fitted.thetas.n_rows != x.n_cols) {
      return missing_standard_errors(x.n_cols);
    }
    arma::colvec const eta = x * fitted.thetas.col(0);
    arma::colvec weights(eta.n_elem);
    if (family == "binomial") {
      for (arma::uword index = 0; index < eta.n_elem; ++index) {
        double const probability = eta(index) >= 0.0
            ? 1.0 / (1.0 + std::exp(-eta(index)))
            : std::exp(eta(index)) / (1.0 + std::exp(eta(index)));
        weights(index) = probability * (1.0 - probability);
      }
    } else {
      weights = arma::exp(eta);
    }
    double const threshold = std::sqrt(std::numeric_limits<double>::epsilon());
    if (!weights.is_finite() || arma::any(weights <= threshold)) {
      return missing_standard_errors(x.n_cols);
    }
    arma::mat const information = x.t() * (x.each_col() % weights);
    double const reciprocal_condition = arma::rcond(information);
    arma::mat inverse;
    if (!std::isfinite(reciprocal_condition) ||
        reciprocal_condition <= threshold ||
        !arma::inv(inverse, information) || !inverse.is_finite() ||
        arma::any(inverse.diag() < 0.0)) {
      return missing_standard_errors(x.n_cols);
    }
    return arma::sqrt(inverse.diag());
  }
  throw std::invalid_argument("fastcpd: Wald intervals are unavailable for " +
                              family);
}

std::vector<ConfidenceInterval> confidence_wald(
    Result const& result, arma::mat const& data,
    ConfidenceOptions const& options) {
  if (result.cp_only || result.thetas.is_empty()) {
    throw std::invalid_argument(
        "fastcpd: Wald intervals require a detailed detector result");
  }
  std::string const family = confidence_native_family(result.family);
  unsigned int const p_response = options.detector_options.p_response;
  if (family == "lm" && p_response > 1) {
    throw std::invalid_argument(
        "fastcpd: Wald intervals are unavailable for multivariate LM");
  }
  arma::mat const analysis_data = confidence_analysis_data(result, data);
  std::vector<arma::uword> const points =
      confidence_change_points(result, data.n_rows);
  std::vector<arma::uword> bounds;
  bounds.reserve(points.size() + 2);
  bounds.push_back(0);
  bounds.insert(bounds.end(), points.begin(), points.end());
  bounds.push_back(data.n_rows);
  if (result.thetas.n_cols != bounds.size() - 1) {
    throw std::invalid_argument(
        "fastcpd: result thetas must have one column per segment");
  }
  double const z = RRandom::standard_normal_quantile(
      1.0 - (1.0 - options.level) / 2.0);
  std::vector<ConfidenceInterval> rows;
  rows.reserve(result.thetas.n_elem);
  for (arma::uword segment = 0; segment < result.thetas.n_cols; ++segment) {
    arma::mat const segment_data =
        analysis_data.rows(bounds[segment], bounds[segment + 1] - 1);
    arma::colvec const standard_errors = theta_standard_errors(
        segment_data, family, result.order, options.detector_options);
    if (standard_errors.n_elem != result.thetas.n_rows) {
      throw std::runtime_error(
          "fastcpd: Wald standard-error dimension differs from theta");
    }
    for (arma::uword parameter = 0; parameter < result.thetas.n_rows;
         ++parameter) {
      ConfidenceInterval row;
      row.parm = "theta";
      row.method = "wald";
      row.segment = static_cast<unsigned int>(segment + 1);
      row.parameter = static_cast<unsigned int>(parameter + 1);
      row.estimate = result.thetas(parameter, segment);
      row.se = standard_errors(parameter);
      row.lower = row.estimate - z * row.se;
      row.upper = row.estimate + z * row.se;
      row.level = options.level;
      rows.push_back(std::move(row));
    }
  }
  return rows;
}

arma::mat segment_bootstrap_data(arma::mat const& data,
                                 std::vector<arma::uword> const& points,
                                 RRandom* random) {
  std::vector<arma::uword> bounds;
  bounds.reserve(points.size() + 2);
  bounds.push_back(0);
  bounds.insert(bounds.end(), points.begin(), points.end());
  bounds.push_back(data.n_rows);
  arma::mat result = data;
  for (std::size_t segment = 0; segment + 1 < bounds.size(); ++segment) {
    arma::uword const start = bounds[segment];
    arma::uword const size = bounds[segment + 1] - start;
    std::vector<arma::uword> const selected =
        random->sample_with_replacement(size, size);
    for (arma::uword offset = 0; offset < size; ++offset) {
      result.row(start + offset) = data.row(start + selected[offset]);
    }
  }
  return result;
}

Result confidence_refit(Result const& result, arma::mat const& data,
                        Options detector_options, RRandom* random) {
  detector_options.family = result.family;
  detector_options.order = result.order;
  detector_options.cp_only = true;
  detector_options.show_progress = false;
  if (result.family == "kcp") {
    return detect_kernel_with_random(data, std::move(detector_options), random);
  }
  return detect(data, std::move(detector_options));
}

double quantile_type1(std::vector<double> values, double probability) {
  std::sort(values.begin(), values.end());
  if (probability <= 0.0) return values.front();
  if (probability >= 1.0) return values.back();
  std::size_t index =
      static_cast<std::size_t>(std::ceil(probability * values.size())) - 1;
  index = std::min(index, values.size() - 1);
  return values[index];
}

std::vector<ConfidenceInterval> confidence_bootstrap(
    Result const& result, arma::mat const& data,
    ConfidenceOptions const& options) {
  if (options.bootstrap != "nonparametric") {
    throw std::invalid_argument(
        "fastcpd: only nonparametric bootstrap is implemented");
  }
  if (options.bootstrap_replicates == 0) {
    throw std::invalid_argument(
        "fastcpd: bootstrap_replicates must be positive");
  }
  std::vector<arma::uword> const points =
      confidence_change_points(result, data.n_rows);
  if (points.empty()) return {};
  RRandom random(kernel_seed(options.seed));
  arma::mat matched(
      options.bootstrap_replicates, points.size(),
      arma::fill::value(std::numeric_limits<double>::quiet_NaN()));
  for (unsigned int replicate = 0; replicate < options.bootstrap_replicates;
       ++replicate) {
    arma::mat const sample = segment_bootstrap_data(data, points, &random);
    std::vector<arma::uword> fitted_points;
    try {
      Result const fitted = confidence_refit(
          result, sample, options.detector_options, &random);
      fitted_points = confidence_change_points(fitted, data.n_rows);
    } catch (...) {
      continue;
    }
    for (std::size_t index = 0; index < points.size(); ++index) {
      arma::uword const left = index == 0
                                  ? 0
                                  : (points[index - 1] + points[index]) / 2;
      arma::uword const right =
          index + 1 == points.size()
              ? data.n_rows
              : static_cast<arma::uword>(std::ceil(
                    (points[index] + points[index + 1]) / 2.0));
      bool found = false;
      arma::uword best = 0;
      arma::uword best_distance = std::numeric_limits<arma::uword>::max();
      for (arma::uword const candidate : fitted_points) {
        if (candidate <= left || candidate > right) continue;
        arma::uword const distance = candidate > points[index]
                                         ? candidate - points[index]
                                         : points[index] - candidate;
        if (!found || distance < best_distance) {
          found = true;
          best = candidate;
          best_distance = distance;
        }
      }
      if (found) matched(replicate, index) = static_cast<double>(best);
    }
  }

  double const alpha = 1.0 - options.level;
  std::vector<ConfidenceInterval> rows;
  rows.reserve(points.size());
  for (std::size_t index = 0; index < points.size(); ++index) {
    std::vector<double> detected;
    for (arma::uword replicate = 0; replicate < matched.n_rows; ++replicate) {
      if (!std::isnan(matched(replicate, index))) {
        detected.push_back(matched(replicate, index));
      }
    }
    ConfidenceInterval row;
    row.parm = "cp";
    row.method = "bootstrap";
    row.bootstrap = options.bootstrap;
    row.index = static_cast<unsigned int>(index + 1);
    row.estimate = static_cast<double>(points[index]);
    row.detection_rate = static_cast<double>(detected.size()) /
                         static_cast<double>(options.bootstrap_replicates);
    row.level = options.level;
    if (!detected.empty()) {
      row.lower = quantile_type1(detected, alpha / 2.0);
      row.upper = quantile_type1(detected, 1.0 - alpha / 2.0);
    }
    rows.push_back(std::move(row));
  }
  return rows;
}

}  // namespace

std::vector<ConfidenceInterval> confint(Result const& result,
                                        arma::mat const& data,
                                        ConfidenceOptions options) {
  require_finite_data(data);
  if (!std::isfinite(options.level) || options.level <= 0.0 ||
      options.level >= 1.0) {
    throw std::invalid_argument("fastcpd: confidence level must be in (0, 1)");
  }
  if (options.parm != "cp" && options.parm != "theta") {
    throw std::invalid_argument("fastcpd: parm must be cp or theta");
  }
  if (options.min_segment_length == 0) {
    throw std::invalid_argument(
        "fastcpd: min_segment_length must be positive");
  }
  if (options.method.empty()) {
    options.method = options.parm == "cp" ? "bootstrap" : "wald";
  }
  if (options.parm == "cp" && options.method == "bootstrap") {
    return confidence_bootstrap(result, data, options);
  }
  if (options.parm == "cp" && options.method == "profile") {
    return confidence_profile(result, data, options);
  }
  if (options.parm == "theta" && options.method == "wald") {
    return confidence_wald(result, data, options);
  }
  throw std::invalid_argument("fastcpd: unavailable confidence method");
}

std::vector<ConfidenceInterval> confint(Result const& result,
                                        arma::colvec const& data,
                                        ConfidenceOptions options) {
  return confint(result, arma::mat(data), std::move(options));
}

}  // namespace fastcpd
