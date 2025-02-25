#ifndef FIT_GLM_H
#define FIT_GLM_H

#include <Rcpp.h>

using namespace Rcpp;

List fastglm(NumericMatrix x,
             SEXP y,
             std::string family,
             Nullable<NumericVector> start = R_NilValue,
             Nullable<NumericVector> weights = R_NilValue,
             Nullable<NumericVector> offset = R_NilValue,
             Nullable<NumericVector> etastart = R_NilValue,
             Nullable<NumericVector> mustart = R_NilValue,
             int method = 0,
             double tol = 1e-8,
             int maxit = 100);

#endif // FIT_GLM_H
