#ifndef REF_TSERIES_H_
#define REF_TSERIES_H_

#include <RcppArmadillo.h>

using ::arma::colvec;
using ::Rcpp::IntegerVector;
using ::Rcpp::List;
using ::Rcpp::Nullable;
using ::Rcpp::NumericVector;
using ::std::string;

List garch(colvec const& x_, colvec const& order_,
           Nullable<string> series = R_NilValue, int maxiter = 200,
           bool trace = false, Nullable<NumericVector> start = R_NilValue,
           string grad = "analytical", Nullable<double> abstol_ = R_NilValue,
           Nullable<double> reltol_ = R_NilValue,
           Nullable<double> xtol_ = R_NilValue,
           Nullable<double> falsetol_ = R_NilValue);

#endif  // REF_TSERIES_H_
