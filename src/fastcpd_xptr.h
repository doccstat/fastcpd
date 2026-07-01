#ifndef FASTCPD_XPTR_H_
#define FASTCPD_XPTR_H_

#include <RcppArmadillo.h>

#include <cstring>
#include <functional>

namespace fastcpd::xptr {

// Canonical raw function-pointer signatures expected from pre-compiled
// (`Rcpp::XPtr`-wrapped) custom cost functions. Users writing compiled cost
// functions define functions matching these signatures, wrap the address in
// an `Rcpp::XPtr<...>` tagged with the matching string constant below, and
// pass the resulting external pointer object through the *existing*
// `cost` / `cost_gradient` / `cost_hessian` arguments of `fastcpd()` --
// exactly where an R closure would otherwise go. See `?fastcpd` for a
// worked example.
using CostPeltFnPtr = double (*)(arma::mat const&);
using CostSenFnPtr = double (*)(arma::mat const&, arma::colvec const&);
using CostGradientFnPtr =
    arma::colvec (*)(arma::mat const&, arma::colvec const&);
using CostHessianFnPtr = arma::mat (*)(arma::mat const&, arma::colvec const&);

constexpr char const* kCostPeltTag = "fastcpd_cost_pelt";
constexpr char const* kCostSenTag = "fastcpd_cost_sen";
constexpr char const* kCostGradientTag = "fastcpd_cost_gradient";
constexpr char const* kCostHessianTag = "fastcpd_cost_hessian";

// Lightweight stand-in for `RcppXPtrUtils::checkXPtrTag` (we avoid that
// dependency entirely): compares the external pointer's R-level tag -- a
// length-one character vector supplied when the pointer was created, e.g.
// `Rcpp::XPtr<CostPeltFnPtr>(new CostPeltFnPtr(&my_cost), true,
// Rcpp::wrap(fastcpd::xptr::kCostPeltTag))` -- against `expected`, throwing a
// clear error on mismatch or absence. This is the package's entire external
// pointer "type system": intentionally minimal, with no external dependency.
inline void CheckXPtrTag(SEXP xptr_sexp, char const* expected) {
  SEXP const tag = R_ExternalPtrTag(xptr_sexp);
  bool const tag_matches = TYPEOF(tag) == STRSXP && Rf_length(tag) == 1 &&
                           std::strcmp(CHAR(STRING_ELT(tag, 0)), expected) == 0;
  if (!tag_matches) {
    Rcpp::stop(
        "Compiled cost function does not carry the expected external "
        "pointer tag \"%s\". Build it with `Rcpp::XPtr<FnPtr>(new "
        "FnPtr(&your_function), true, Rcpp::wrap(\"%s\"))` -- see `?fastcpd` "
        "for the exact function signature each argument expects.",
        expected, expected);
  }
}

inline std::function<double(arma::mat const&)> WrapCostPelt(SEXP xptr_sexp) {
  CheckXPtrTag(xptr_sexp, kCostPeltTag);
  CostPeltFnPtr const fn = *Rcpp::XPtr<CostPeltFnPtr>(xptr_sexp);
  return [fn](arma::mat const& data) -> double { return fn(data); };
}

inline std::function<double(arma::mat const&, arma::colvec const&)> WrapCostSen(
    SEXP xptr_sexp) {
  CheckXPtrTag(xptr_sexp, kCostSenTag);
  CostSenFnPtr const fn = *Rcpp::XPtr<CostSenFnPtr>(xptr_sexp);
  return [fn](arma::mat const& data, arma::colvec const& theta) -> double {
    return fn(data, theta);
  };
}

inline std::function<arma::colvec(arma::mat const&, arma::colvec const&)> WrapCostGradient(
    SEXP xptr_sexp) {
  CheckXPtrTag(xptr_sexp, kCostGradientTag);
  CostGradientFnPtr const fn = *Rcpp::XPtr<CostGradientFnPtr>(xptr_sexp);
  return [fn](arma::mat const& data, arma::colvec const& theta) -> arma::colvec {
    return fn(data, theta);
  };
}

inline std::function<arma::mat(arma::mat const&, arma::colvec const&)> WrapCostHessian(
    SEXP xptr_sexp) {
  CheckXPtrTag(xptr_sexp, kCostHessianTag);
  CostHessianFnPtr const fn = *Rcpp::XPtr<CostHessianFnPtr>(xptr_sexp);
  return [fn](arma::mat const& data, arma::colvec const& theta) -> arma::mat {
    return fn(data, theta);
  };
}

}  // namespace fastcpd::xptr

#endif  // FASTCPD_XPTR_H_
