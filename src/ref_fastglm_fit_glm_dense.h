#ifndef FIT_GLM_DENSE_H
#define FIT_GLM_DENSE_H

#include <Rcpp.h>

using namespace Rcpp;

List fit_glm(Rcpp::NumericMatrix Xs,
             Rcpp::NumericVector ys,
             Rcpp::NumericVector weightss,
             Rcpp::NumericVector offsets,
             Rcpp::NumericVector starts,
             Rcpp::NumericVector etas,
             int type,
             double tol,
             int maxit,
             std::string family);

#endif // FIT_GLM_DENSE_H
