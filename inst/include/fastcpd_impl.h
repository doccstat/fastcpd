#ifndef FASTCPD_IMPL_H_
#define FASTCPD_IMPL_H_

#include "RcppArmadillo.h"

using ::arma::colvec;
using ::arma::mat;
using ::Rcpp::Function;
using ::Rcpp::List;
using ::Rcpp::Nullable;
using ::std::function;
using ::std::string;

List fastcpd_impl(
    mat data,
    double beta,
    const string cost_adjustment,
    const unsigned int d,
    const int segment_count,
    const double trim,
    const double momentum_coef,
    Nullable<Function> k,
    const string family,
    const double epsilon,
    const int p,
    const colvec order,
    Nullable<Function> cost,
    Nullable<Function> cost_gradient,
    Nullable<Function> cost_hessian,
    const bool cp_only,
    const double vanilla_percentage,
    const bool warm_start,
    colvec lower,
    colvec upper,
    colvec line_search,
    const mat variance_estimate,
    const unsigned int p_response,
    const double pruning_coef,
    const string r_clock,
    const bool r_progress
);

#endif  // FASTCPD_IMPL_H_
