#ifndef FIT_GLM_DENSE_H
#define FIT_GLM_DENSE_H

#include <RcppArmadillo.h>

using ::Rcpp::List;
using ::Rcpp::NumericMatrix;
using ::Rcpp::NumericVector;
using ::std::string;

List fit_glm(NumericMatrix Xs, NumericVector ys, NumericVector weightss,
             NumericVector offsets, NumericVector starts, NumericVector etas,
             int type, double tol, int maxit, string family);

#endif  // FIT_GLM_DENSE_H
