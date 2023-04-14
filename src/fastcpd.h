#ifndef OCPD_H
#define OCPD_H

#include <RcppArmadillo.h>
#include <testthat.h>

arma::colvec cost_update_gradient(arma::mat, arma::colvec, std::string);
arma::mat cost_update_hessian(arma::mat, arma::colvec, std::string, double);
Rcpp::List cost_update_c(
    const arma::cube,
    arma::mat,
    arma::mat,
    arma::cube,
    const int,
    const int,
    const int,
    const std::string,
    arma::colvec,
    const double,
    const double,
    const double,
    const double,
    const double,
    const double,
    Rcpp::Function,
    Rcpp::Function
);
double negative_log_likelihood_c(arma::mat, arma::colvec, std::string, double, bool);

#endif
