#ifndef NO_RCPP
#define NO_RCPP
#endif

#include <fastcpd/fastcpd.h>

#include "fastcpd_template.h"
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
         family == "mgaussian" || family == "garch";
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
      (family == "mgaussian" || family == "lasso" || family == "garch")) {
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

arma::mat nearest_positive_definite(arma::mat const& matrix) {
  if (matrix.n_rows != matrix.n_cols) {
    throw std::invalid_argument(
        "fastcpd: variance_estimate must be a square matrix");
  }
  if (matrix.n_rows == 0) return matrix;
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  arma::eig_sym(eigenvalues, eigenvectors, arma::symmatu(matrix));
  double const floor = std::max(1e-12, std::numeric_limits<double>::epsilon());
  eigenvalues.transform([floor](double value) {
    return value > floor ? value : floor;
  });
  return eigenvectors * arma::diagmat(eigenvalues) * eigenvectors.t();
}

arma::mat estimate_mean_variance(arma::mat const& data) {
  if (data.n_rows < 2) return arma::eye(data.n_cols, data.n_cols);
  arma::mat const diffs =
      data.rows(1, data.n_rows - 1) - data.rows(0, data.n_rows - 2);
  return nearest_positive_definite(diffs.t() * diffs /
                                   (2.0 * static_cast<double>(diffs.n_rows)));
}

arma::mat estimate_mgaussian_variance(arma::mat const& data,
                                      unsigned int p_response) {
  unsigned int const q = p_response > 0 ? p_response : data.n_cols;
  if (q == 0 || q > data.n_cols) {
    throw std::invalid_argument(
        "fastcpd: p_response must be between 1 and data.n_cols");
  }
  arma::mat const y = data.cols(0, q - 1);
  arma::mat residuals;
  if (q < data.n_cols) {
    arma::mat const x = data.cols(q, data.n_cols - 1);
    residuals = y - x * (arma::pinv(x) * y);
  } else {
    residuals = y.each_row() - arma::mean(y, 0);
  }
  if (residuals.n_rows < 2) return arma::eye(q, q);
  arma::mat const diffs =
      residuals.rows(1, residuals.n_rows - 1) -
      residuals.rows(0, residuals.n_rows - 2);
  return nearest_positive_definite(diffs.t() * diffs /
                                   (2.0 * static_cast<double>(diffs.n_rows)));
}

arma::mat resolve_variance_estimate(Options const& options,
                                    arma::mat const& data,
                                    std::string const& family) {
  if (!options.variance_estimate.is_empty()) {
    arma::mat value = nearest_positive_definite(options.variance_estimate);
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
    return value;
  }
  if (family == "mean") return estimate_mean_variance(data);
  if (family == "mgaussian") {
    return estimate_mgaussian_variance(data, options.p_response);
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
  return Result{std::get<0>(value), std::get<1>(value), std::get<2>(value),
                std::get<3>(value), std::get<4>(value)};
}

}  // namespace detail

Result detect(arma::mat const& data, Options options) {
  if (data.n_rows == 0 || data.n_cols == 0) {
    throw std::invalid_argument("fastcpd: data must be a non-empty matrix");
  }
  detail::validate_cost_adjustment(options.cost_adjustment);
  std::string const family =
      detail::normalize_family(options.family, options.order);
  detail::validate_custom_options(options, family);

  arma::colvec const order = options.order;
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
  arma::mat const variance_estimate =
      detail::resolve_variance_estimate(options, data, family);

  double const beta =
      detail::compute_beta(options.beta, options.beta_criterion,
                           static_cast<int>(data.n_rows), p);
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
