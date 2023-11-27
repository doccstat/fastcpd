#ifndef FASTCPD_TYPES_H_
#define FASTCPD_TYPES_H_

#include "RcppArmadillo.h"

using ::arma::accu;
using ::arma::as_scalar;
using ::arma::clamp;
using ::arma::eye;
using ::arma::index_min;
using ::arma::join_cols;
using ::arma::linspace;
using ::arma::log_det_sympd;
using ::arma::max;
using ::arma::min;
using ::arma::norm;
using ::arma::ones;
using ::arma::sign;
using ::arma::solve;
using ::arma::zeros;
using ::Rcpp::as;
using ::Rcpp::Named;
using ::Rcpp::stop;

using ::arma::colvec;
using ::arma::cube;
using ::arma::mat;
using ::arma::rowvec;
using ::arma::ucolvec;
using ::arma::uvec;
using ::arma::vec;
using ::Rcpp::Environment;
using ::Rcpp::Function;
using ::Rcpp::InternalFunction;
using ::Rcpp::List;
using ::Rcpp::Nullable;
using ::Rcpp::NumericVector;
using ::Rcpp::S4;
using ::std::function;
using ::std::string;

#endif  // FASTCPD_TYPES_H_
