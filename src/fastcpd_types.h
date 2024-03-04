#ifndef FASTCPD_TYPES_H_
#define FASTCPD_TYPES_H_

#include "RcppArmadillo.h"

#ifdef DEBUG
#  define DEBUG_RCOUT(x) Rcout << #x << ": " << x << std::endl
#else
#  define DEBUG_RCOUT(x) do {} while (0)
#endif

using ::arma::abs;
using ::arma::accu;
using ::arma::approx_equal;
using ::arma::as_scalar;
using ::arma::clamp;
using ::arma::cov;
using ::arma::diff;
using ::arma::dot;
using ::arma::eye;
using ::arma::floor;
using ::arma::index_min;
using ::arma::join_cols;
using ::arma::join_rows;
using ::arma::join_slices;
using ::arma::linspace;
using ::arma::log_det_sympd;
using ::arma::max;
using ::arma::mean;
using ::arma::min;
using ::arma::norm;
using ::arma::ones;
using ::arma::reverse;
using ::arma::sign;
using ::arma::solve;
using ::arma::sort;
using ::arma::square;
using ::arma::trace;
using ::arma::unique;
using ::arma::zeros;
using ::Rcpp::as;
using ::Rcpp::Named;
using ::Rcpp::stop;
using ::Rcpp::wrap;

using ::arma::colvec;
using ::arma::cube;
using ::arma::mat;
using ::arma::rowvec;
using ::arma::span;
using ::arma::ucolvec;
using ::arma::uvec;
using ::arma::vec;
using ::Rcpp::checkUserInterrupt;
using ::Rcpp::Environment;
using ::Rcpp::Function;
using ::Rcpp::InternalFunction;
using ::Rcpp::List;
using ::Rcpp::Nullable;
using ::Rcpp::NumericVector;
using ::Rcpp::Rcout;
using ::Rcpp::S4;
using ::std::function;
using ::std::string;

#endif  // FASTCPD_TYPES_H_
