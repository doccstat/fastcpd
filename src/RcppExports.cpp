// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "fastcpd_types.h"
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fastcpd_impl
List fastcpd_impl(const mat data, const double beta, const string cost_adjustment, const unsigned int d, const int segment_count, const double trim, const double momentum_coef, const Nullable<Function> multiple_epochs_function, const string family, const double epsilon, const int p, const colvec order, const Nullable<Function> cost, const Nullable<Function> cost_gradient, const Nullable<Function> cost_hessian, const bool cp_only, const double vanilla_percentage, const bool warm_start, const colvec lower, const colvec upper, const colvec line_search, const mat variance_estimate, const unsigned int p_response, const double pruning_coef, const string r_clock, const bool r_progress);
RcppExport SEXP _fastcpd_fastcpd_impl(SEXP dataSEXP, SEXP betaSEXP, SEXP cost_adjustmentSEXP, SEXP dSEXP, SEXP segment_countSEXP, SEXP trimSEXP, SEXP momentum_coefSEXP, SEXP multiple_epochs_functionSEXP, SEXP familySEXP, SEXP epsilonSEXP, SEXP pSEXP, SEXP orderSEXP, SEXP costSEXP, SEXP cost_gradientSEXP, SEXP cost_hessianSEXP, SEXP cp_onlySEXP, SEXP vanilla_percentageSEXP, SEXP warm_startSEXP, SEXP lowerSEXP, SEXP upperSEXP, SEXP line_searchSEXP, SEXP variance_estimateSEXP, SEXP p_responseSEXP, SEXP pruning_coefSEXP, SEXP r_clockSEXP, SEXP r_progressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const string >::type cost_adjustment(cost_adjustmentSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type d(dSEXP);
    Rcpp::traits::input_parameter< const int >::type segment_count(segment_countSEXP);
    Rcpp::traits::input_parameter< const double >::type trim(trimSEXP);
    Rcpp::traits::input_parameter< const double >::type momentum_coef(momentum_coefSEXP);
    Rcpp::traits::input_parameter< const Nullable<Function> >::type multiple_epochs_function(multiple_epochs_functionSEXP);
    Rcpp::traits::input_parameter< const string >::type family(familySEXP);
    Rcpp::traits::input_parameter< const double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const colvec >::type order(orderSEXP);
    Rcpp::traits::input_parameter< const Nullable<Function> >::type cost(costSEXP);
    Rcpp::traits::input_parameter< const Nullable<Function> >::type cost_gradient(cost_gradientSEXP);
    Rcpp::traits::input_parameter< const Nullable<Function> >::type cost_hessian(cost_hessianSEXP);
    Rcpp::traits::input_parameter< const bool >::type cp_only(cp_onlySEXP);
    Rcpp::traits::input_parameter< const double >::type vanilla_percentage(vanilla_percentageSEXP);
    Rcpp::traits::input_parameter< const bool >::type warm_start(warm_startSEXP);
    Rcpp::traits::input_parameter< const colvec >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< const colvec >::type upper(upperSEXP);
    Rcpp::traits::input_parameter< const colvec >::type line_search(line_searchSEXP);
    Rcpp::traits::input_parameter< const mat >::type variance_estimate(variance_estimateSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type p_response(p_responseSEXP);
    Rcpp::traits::input_parameter< const double >::type pruning_coef(pruning_coefSEXP);
    Rcpp::traits::input_parameter< const string >::type r_clock(r_clockSEXP);
    Rcpp::traits::input_parameter< const bool >::type r_progress(r_progressSEXP);
    rcpp_result_gen = Rcpp::wrap(fastcpd_impl(data, beta, cost_adjustment, d, segment_count, trim, momentum_coef, multiple_epochs_function, family, epsilon, p, order, cost, cost_gradient, cost_hessian, cp_only, vanilla_percentage, warm_start, lower, upper, line_search, variance_estimate, p_response, pruning_coef, r_clock, r_progress));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP run_testthat_tests(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_fastcpd_fastcpd_impl", (DL_FUNC) &_fastcpd_fastcpd_impl, 26},
    {"run_testthat_tests", (DL_FUNC) &run_testthat_tests, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_fastcpd(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
