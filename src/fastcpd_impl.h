#ifndef FASTCPD_IMPL_H_
#define FASTCPD_IMPL_H_

#include "RcppArmadillo.h"

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
    bool const r_progress);

#endif  // FASTCPD_IMPL_H_
