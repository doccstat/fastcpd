#ifndef FASTCPD_IMPL_H_
#define FASTCPD_IMPL_H_

#include "RcppArmadillo.h"

Rcpp::List fastcpd_impl(
    const arma::mat& data, const double beta,
    const std::string& cost_adjustment, const int segment_count,
    const double trim, const double momentum_coef,
    const Rcpp::Nullable<Rcpp::Function>& multiple_epochs_function,
    const std::string& family, const double epsilon, const int p,
    const arma::colvec& order, const Rcpp::Nullable<Rcpp::Function>& cost,
    const Rcpp::Nullable<Rcpp::Function>& cost_gradient,
    const Rcpp::Nullable<Rcpp::Function>& cost_hessian, const bool cp_only,
    const double vanilla_percentage, const bool warm_start,
    const arma::colvec& lower, const arma::colvec& upper,
    const arma::colvec& line_search, const arma::mat& variance_estimate,
    const unsigned int p_response, const double pruning_coef,
    const bool r_progress);

#endif  // FASTCPD_IMPL_H_
