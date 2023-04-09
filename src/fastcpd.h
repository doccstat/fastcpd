#ifndef OCPD_H
#define OCPD_H

#include <RcppArmadillo.h>
#include <testthat.h>

arma::colvec cost_update_gradient(arma::rowvec, arma::colvec, std::string);
arma::mat cost_update_hessian(arma::rowvec, arma::colvec, std::string, double);

#endif
