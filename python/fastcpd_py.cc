// NO_RCPP must be defined before any fastcpd header is included so that
// fastcpd_template.h uses <armadillo> instead of <RcppArmadillo.h>.
#define NO_RCPP

#include "fastcpd_py.h"
#include "fastcpd_template.h"
#include "families/lasso.h"
#include "families/mean.h"
#include "families/meanvariance.h"
#include "families/mgaussian.h"
#include "families/variance.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace fastcpd {

// ---------------------------------------------------------------------------
// Helper implementations
// ---------------------------------------------------------------------------

int fastcpd_compute_p(std::string const& family, int n_cols,
                      int p_response, int p_explicit) {
  if (p_explicit > 0) return p_explicit;
  if (family == "mean")         return n_cols;
  if (family == "variance")     return n_cols * n_cols;
  if (family == "meanvariance") return n_cols + n_cols * n_cols;
  if (family == "mgaussian") {
    int q = (p_response > 0) ? p_response : n_cols;
    int k = n_cols - q;
    return (k > 0) ? q * k : q;
  }
  if (family == "lasso") return n_cols - 1;
  return n_cols;  // default fallback
}

double fastcpd_compute_beta(double beta_val, std::string const& beta_str,
                            int n, int p) {
  if (beta_str.empty()) return beta_val;
  double const ln    = std::log(static_cast<double>(n));
  double const log2n = ln / std::log(2.0);
  if (beta_str == "BIC")  return (p + 1) * ln    / 2.0;
  if (beta_str == "MBIC") return (p + 2) * ln    / 2.0;
  if (beta_str == "MDL")  return (p + 2) * log2n / 2.0;
  throw std::runtime_error(
      "fastcpd: unknown beta criterion '" + beta_str +
      "'. Use \"BIC\", \"MBIC\", \"MDL\", or supply a numeric value.");
}

double fastcpd_compute_pruning_coef(double pruning_coef,
                                    bool pruning_coef_explicitly_set,
                                    std::string const& cost_adjustment,
                                    std::string const& family, int p) {
  if (pruning_coef_explicitly_set) return pruning_coef;
  // No pruning for families that are sensitive to over-pruning.
  if (family == "mgaussian" || family == "lasso")
    pruning_coef = -std::numeric_limits<double>::infinity();
  // Add penalty adjustment.  Adding finite to -Inf stays -Inf (no pruning).
  if (cost_adjustment == "MBIC")
    pruning_coef += p * std::log(2.0);
  else if (cost_adjustment == "MDL")
    pruning_coef += p * std::log2(2.0);
  return pruning_coef;
}

// ---------------------------------------------------------------------------
// Internal dispatch
// ---------------------------------------------------------------------------

namespace {

using namespace fastcpd::classes;
using Result =
    std::tuple<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat>;

// Convenience alias for the full argument list forwarded to Fastcpd<>.
#define FWD_ARGS \
    beta, nullptr, nullptr, nullptr, nullptr, cp_only, data, epsilon, family, \
    nullptr, line_search, lower, momentum_coef, order, p, p_response,        \
    pruning_coef, segment_count, trim, upper, vanilla_percentage,             \
    variance_estimate, warm_start

// Instantiate one Fastcpd<> variant, run it, and return the result.
// kNDims: -1 = general path, 1 = scalar fast path.
#define MAKE_AND_RUN(Policy, kCost, kNDims_val)                            \
  {                                                                        \
    fastcpd::classes::Fastcpd<fastcpd::families::Policy, false, true,     \
                              fastcpd::classes::CostAdjustment::kCost,    \
                              false, kNDims_val>                           \
        inst(FWD_ARGS);                                                    \
    return inst.Run();                                                     \
  }

// Dispatch on cost_adjustment for a family that has a 1-D fast path.
#define DISPATCH_PELT_1D(Policy)                                           \
  if (data.n_cols == 1) {                                                  \
    if (cost_adjustment == "MBIC") MAKE_AND_RUN(Policy, kMBIC, 1);        \
    if (cost_adjustment == "MDL")  MAKE_AND_RUN(Policy, kMDL,  1);        \
    MAKE_AND_RUN(Policy, kBIC, 1);                                         \
  } else {                                                                 \
    if (cost_adjustment == "MBIC") MAKE_AND_RUN(Policy, kMBIC, -1);       \
    if (cost_adjustment == "MDL")  MAKE_AND_RUN(Policy, kMDL,  -1);       \
    MAKE_AND_RUN(Policy, kBIC, -1);                                        \
  }

// Dispatch on cost_adjustment for a family without a 1-D fast path.
#define DISPATCH_FAMILY(Policy)                                            \
  if (cost_adjustment == "MBIC") MAKE_AND_RUN(Policy, kMBIC, -1);        \
  if (cost_adjustment == "MDL")  MAKE_AND_RUN(Policy, kMDL,  -1);        \
  MAKE_AND_RUN(Policy, kBIC, -1);

// --- forwarding parameters (match Fastcpd constructor under NO_RCPP) ---
Result dispatch_impl(
    double beta, bool cp_only, arma::mat const& data, double epsilon,
    std::string const& family, arma::colvec const& line_search,
    arma::colvec const& lower, double momentum_coef, arma::colvec const& order,
    int p, unsigned int p_response, double pruning_coef, int segment_count,
    double trim, arma::colvec const& upper, double vanilla_percentage,
    arma::mat const& variance_estimate, bool warm_start,
    std::string const& cost_adjustment) {

  if (family == "mean")        { DISPATCH_PELT_1D(MeanFamily); }
  if (family == "variance")    { DISPATCH_PELT_1D(VarianceFamily); }
  if (family == "meanvariance"){ DISPATCH_PELT_1D(MeanvarianceFamily); }
  if (family == "mgaussian")   { DISPATCH_FAMILY(MgaussianFamily); }
  if (family == "lasso")       { DISPATCH_FAMILY(LassoFamily); }

  throw std::runtime_error(
      "fastcpd Python: unsupported family '" + family +
      "'. Supported: mean, variance, meanvariance, mgaussian, lasso.");
}

#undef FWD_ARGS
#undef MAKE_AND_RUN
#undef DISPATCH_PELT_1D
#undef DISPATCH_FAMILY

}  // namespace

std::tuple<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat>
fastcpd_py_dispatch(
    double beta, std::string const& cost_adjustment, bool cp_only,
    arma::mat const& data, double epsilon, std::string const& family,
    arma::colvec const& line_search, arma::colvec const& lower,
    double momentum_coef, arma::colvec const& order, int p,
    unsigned int p_response, double pruning_coef, int segment_count,
    double trim, arma::colvec const& upper, double vanilla_percentage,
    arma::mat const& variance_estimate, bool warm_start) {
  return dispatch_impl(beta, cp_only, data, epsilon, family, line_search,
                       lower, momentum_coef, order, p, p_response,
                       pruning_coef, segment_count, trim, upper,
                       vanilla_percentage, variance_estimate, warm_start,
                       cost_adjustment);
}

}  // namespace fastcpd
