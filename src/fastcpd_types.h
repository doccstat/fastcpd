#ifndef FASTCPD_TYPES_H_
#define FASTCPD_TYPES_H_

#define ARMA_64BIT_WORD

#include "RcppArmadillo.h"

#define ERROR(msg) \
  Rcout << "error: " << __FILE__ << ":" << __LINE__ << ": " << msg << std::endl
#ifdef DEBUG
#define DEBUG_RCOUT(x) Rcout << #x << ": " << x << std::endl
#define INFO(msg) \
  Rcout << "info: " << __FILE__ << ":" << __LINE__ << ": " << msg << std::endl
#else
#define DEBUG_RCOUT(x) do {} while (0)
#define INFO(msg) do {} while (0)
#endif

using ::arma::abs;
using ::arma::accu;
using ::arma::approx_equal;
using ::arma::as_scalar;
using ::arma::clamp;
using ::arma::cov;
using ::arma::diff;
using ::arma::dot;
using ::arma::find;
using ::arma::eye;
using ::arma::floor;
using ::arma::index_max;
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
using ::Rcpp::checkUserInterrupt;
using ::Rcpp::Named;
using ::Rcpp::stop;
using ::Rcpp::wrap;

using ::arma::colvec;
using ::arma::conv_to;
using ::arma::cube;
using ::arma::mat;
using ::arma::rowvec;
using ::arma::span;
using ::arma::ucolvec;
using ::arma::uvec;
using ::arma::vec;
using ::Rcpp::Environment;
using ::Rcpp::Function;
using ::Rcpp::InternalFunction;
using ::Rcpp::List;
using ::Rcpp::Nullable;
using ::Rcpp::NumericMatrix;
using ::Rcpp::NumericVector;
using ::Rcpp::S4;
using ::std::function;
using ::std::string;

using ::Rcpp::Rcout;

#endif  // FASTCPD_TYPES_H_
