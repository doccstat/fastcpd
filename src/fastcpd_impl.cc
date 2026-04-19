#include "fastcpd_impl.h"

// Signature shared by every generated fastcpd_*_impl function.
// r_progress, cost_adjustment, vanilla_percentage==1, and line_search!=c(1)
// have all been consumed by dispatch below.
#define IMPL_PARAMS \
    arma::mat const& data, double const beta, \
    int const segment_count, \
    double const trim, double const momentum_coef, \
    Rcpp::Nullable<Rcpp::Function> const& multiple_epochs_function, \
    std::string const& family, double const epsilon, int const p, \
    arma::colvec const& order, Rcpp::Nullable<Rcpp::Function> const& cost_pelt, \
    Rcpp::Nullable<Rcpp::Function> const& cost_sen, \
    Rcpp::Nullable<Rcpp::Function> const& cost_gradient, \
    Rcpp::Nullable<Rcpp::Function> const& cost_hessian, bool const cp_only, \
    double const vanilla_percentage, bool const warm_start, \
    arma::colvec const& lower, arma::colvec const& upper, \
    arma::colvec const& line_search, arma::mat const& variance_estimate, \
    unsigned int const p_response, double const pruning_coef

// PELT families: 3 cost × 2 progress = 6 variants each.
#define DECLARE_PELT(name) \
    Rcpp::List fastcpd_##name##_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mbic_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mdl_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_prog_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mbic_prog_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mdl_prog_impl(IMPL_PARAMS);

// SEGD families: 3 cost × 2 vanilla × 2 line_search × 2 progress = 24 variants each.
#define DECLARE_SEGD(name) \
    Rcpp::List fastcpd_##name##_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_ls_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_van_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_van_ls_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mbic_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mbic_ls_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mbic_van_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mbic_van_ls_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mdl_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mdl_ls_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mdl_van_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mdl_van_ls_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_prog_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_ls_prog_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_van_prog_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_van_ls_prog_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mbic_prog_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mbic_ls_prog_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mbic_van_prog_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mbic_van_ls_prog_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mdl_prog_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mdl_ls_prog_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mdl_van_prog_impl(IMPL_PARAMS); \
    Rcpp::List fastcpd_##name##_mdl_van_ls_prog_impl(IMPL_PARAMS);

DECLARE_PELT(mean)
DECLARE_PELT(mean_1d)
DECLARE_PELT(mgaussian)
DECLARE_PELT(variance)
DECLARE_PELT(variance_1d)
DECLARE_PELT(meanvariance)
DECLARE_PELT(meanvariance_1d)
DECLARE_PELT(garch)
DECLARE_SEGD(arma)
DECLARE_SEGD(binomial)
DECLARE_SEGD(custom)
DECLARE_SEGD(gaussian)
DECLARE_SEGD(lasso)
DECLARE_SEGD(ma)
DECLARE_SEGD(poisson)

#undef DECLARE_PELT
#undef DECLARE_SEGD

#define ARGS \
    data, beta, segment_count, trim, momentum_coef, \
    multiple_epochs_function, family, epsilon, p, order, cost_pelt, cost_sen, \
    cost_gradient, cost_hessian, cp_only, vanilla_percentage, warm_start, \
    lower, upper, line_search, variance_estimate, p_response, pruning_coef

// Dispatch on r_progress for a fully-qualified base name.
#define CALL_PROG(base) \
    return r_progress ? base##_prog_impl(ARGS) : base##_impl(ARGS)

// Dispatch on cost_adjustment for PELT families.
#define DISPATCH_PELT(name) \
    if (cost_adjustment == "MBIC") { CALL_PROG(fastcpd_##name##_mbic); } \
    else if (cost_adjustment == "MDL") { CALL_PROG(fastcpd_##name##_mdl); } \
    else { CALL_PROG(fastcpd_##name); }

// Dispatch on cost_adjustment + vanilla + line_search for SEGD families.
// Suffix order: [_mbic|_mdl][_van][_ls]
#define DISPATCH_SEGD_LS(base_with_van) \
    if (ls) { CALL_PROG(base_with_van##_ls); } \
    else     { CALL_PROG(base_with_van); }

#define DISPATCH_SEGD_VAN(base_with_cost) \
    if (van) { DISPATCH_SEGD_LS(base_with_cost##_van); } \
    else      { DISPATCH_SEGD_LS(base_with_cost); }

#define DISPATCH_SEGD(name) \
    bool const van = (vanilla_percentage == 1); \
    bool const ls  = (line_search.n_elem > 1 || line_search(0) != 1.0); \
    if (cost_adjustment == "MBIC") { DISPATCH_SEGD_VAN(fastcpd_##name##_mbic); } \
    else if (cost_adjustment == "MDL") { DISPATCH_SEGD_VAN(fastcpd_##name##_mdl); } \
    else { DISPATCH_SEGD_VAN(fastcpd_##name); }

// [[Rcpp::export]]
Rcpp::List fastcpd_impl(
    arma::mat const& data, double const beta,
    std::string const& cost_adjustment, int const segment_count,
    double const trim, double const momentum_coef,
    Rcpp::Nullable<Rcpp::Function> const& multiple_epochs_function,
    std::string const& family, double const epsilon, int const p,
    arma::colvec const& order, Rcpp::Nullable<Rcpp::Function> const& cost_pelt,
    Rcpp::Nullable<Rcpp::Function> const& cost_sen,
    Rcpp::Nullable<Rcpp::Function> const& cost_gradient,
    Rcpp::Nullable<Rcpp::Function> const& cost_hessian, bool const cp_only,
    double const vanilla_percentage, bool const warm_start,
    arma::colvec const& lower, arma::colvec const& upper,
    arma::colvec const& line_search, arma::mat const& variance_estimate,
    unsigned int const p_response, double const pruning_coef,
    bool const r_progress) {

  if (family == "mean") {
    if (data.n_cols == 1) { DISPATCH_PELT(mean_1d); }
    else                  { DISPATCH_PELT(mean); }
  }
  if (family == "mgaussian")                     { DISPATCH_PELT(mgaussian); }
  if (family == "variance") {
    if (data.n_cols == 1) { DISPATCH_PELT(variance_1d); }
    else                  { DISPATCH_PELT(variance); }
  }
  if (family == "meanvariance") {
    if (data.n_cols == 1) { DISPATCH_PELT(meanvariance_1d); }
    else                  { DISPATCH_PELT(meanvariance); }
  }
  if (family == "garch")                         { DISPATCH_PELT(garch); }
  if (family == "arma" || family == "arima")     { DISPATCH_SEGD(arma); }
  if (family == "binomial")                      { DISPATCH_SEGD(binomial); }
  if (family == "gaussian" || family == "lm")    { DISPATCH_SEGD(gaussian); }
  if (family == "lasso")                         { DISPATCH_SEGD(lasso); }
  if (family == "ma")                            { DISPATCH_SEGD(ma); }
  if (family == "poisson")                       { DISPATCH_SEGD(poisson); }
  DISPATCH_SEGD(custom);
}

#undef IMPL_PARAMS
#undef ARGS
#undef CALL_PROG
#undef DISPATCH_PELT
#undef DISPATCH_SEGD_LS
#undef DISPATCH_SEGD_VAN
#undef DISPATCH_SEGD
