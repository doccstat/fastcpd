#ifndef OCPD_H
#define OCPD_H

#include <RcppArmadillo.h>
#include <testthat.h>

arma::colvec cost_update_gradient(arma::mat, arma::colvec, std::string);
arma::mat cost_update_hessian(arma::mat, arma::colvec, std::string, double);
Rcpp::List cost_update(const arma::cube, arma::mat, arma::mat, arma::cube, const int, const int, const int, const std::string, arma::colvec, const double, const double, const double, const double, const double, const double, Rcpp::Function, Rcpp::Function);
Rcpp::List negative_log_likelihood_c(arma::mat, Rcpp::Nullable<arma::colvec>, std::string, double, bool);

#endif
