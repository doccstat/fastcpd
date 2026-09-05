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
#include <cctype>
#include <cmath>
#include <limits>
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
                                              unsigned int d) {
  if (d == 0 || d > data.n_cols) {
    throw std::invalid_argument(
        "fastcpd: p_response must be between 1 and data.n_cols");
  }
  unsigned int const predictor_count =
      static_cast<unsigned int>(data.n_cols - d);
  // R's default block size is predictor_count + 1.  With no predictors the
  // regression solve is 0×0 and variance.lm() yields no usable estimates;
  // callers handle this degenerate direct-mgaussian case below.
  unsigned int const block_size = predictor_count + 1;
  if (predictor_count == 0 || data.n_rows <= block_size) {
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
    double const threshold = q75 + std::numeric_limits<double>::infinity() *
                                       (q75 - q25);
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

Result detect(arma::mat const& data, Options options) {
  if (data.n_rows == 0 || data.n_cols == 0) {
    throw std::invalid_argument("fastcpd: data must be a non-empty matrix");
  }
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

  return detail::make_result(detail::dispatch(
      beta, options.cost_adjustment, options, data, family, order, p,
      p_response, pruning_coef, options.segment_count, vanilla_percentage,
      variance_estimate, lower, upper, line_search));
}

Result detect(arma::colvec const& data, Options options) {
  return detect(arma::mat(data), std::move(options));
}

Result mean(arma::mat const& data, Options options) {
  options.family = "mean";
  return detect(data, std::move(options));
}

Result variance(arma::mat const& data, Options options) {
  options.family = "variance";
  return detect(data, std::move(options));
}

Result meanvariance(arma::mat const& data, Options options) {
  options.family = "meanvariance";
  return detect(data, std::move(options));
}

Result exponential(arma::mat const& data, Options options) {
  options.family = "exponential";
  return detect(data, std::move(options));
}

Result gaussian(arma::mat const& data, Options options) {
  options.family = "gaussian";
  return detect(data, std::move(options));
}

}  // namespace fastcpd
