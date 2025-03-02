#ifndef FIT_GLM_H
#define FIT_GLM_H

#include <RcppArmadillo.h>
#include <RcppEigen.h>

using ::arma::colvec;
using ::arma::mat;
using ::Eigen::Map;
using ::Eigen::VectorXd;
using ::Rcpp::List;
using ::Rcpp::Nullable;
using ::Rcpp::NumericMatrix;
using ::Rcpp::NumericVector;
using ::std::string;

bool valideta_gaussian(const VectorXd &eta);

bool valideta_binomial(const VectorXd &eta);

bool valideta_poisson(const VectorXd &eta);

bool validmu_gaussian(const VectorXd &mu);

bool validmu_binomial(const VectorXd &mu);

bool validmu_poisson(const VectorXd &mu);

// Gaussian deviance residuals: wt * ((y - mu)^2)
NumericVector dev_resids_gaussian(const Map<VectorXd> &y, const VectorXd &mu,
                                  const Map<VectorXd> &wt);

// Binomial deviance residuals: simply call the exported C function.
// Note: binomial_dev_resids_cpp is assumed to wrap the SEXP
// version defined elsewhere.
NumericVector dev_resids_binomial(const Map<VectorXd> &y, const VectorXd &mu,
                                  const Map<VectorXd> &wt);

// Poisson deviance residuals:
//   r = mu * wt,
//   for indices where y > 0, set r = wt * (y * log(y/mu) - (y - mu))
//   and return 2 * r.
NumericVector dev_resids_poisson(const Map<VectorXd> &y, const VectorXd &mu,
                                 const Map<VectorXd> &wt);

// Gaussian variance: always 1 for each element.
NumericVector var_gaussian(const VectorXd &mu);

// Binomial variance: mu * (1 - mu)
NumericVector var_binomial(const VectorXd &mu);

// Poisson variance: just mu
NumericVector var_poisson(const VectorXd &mu);

// Gaussian link inverse: identity (eta).
NumericVector linkinv_gaussian(const VectorXd &eta);

// Binomial link inverse: delegate to logit_linkinv_cpp.
NumericVector linkinv_binomial(const VectorXd &eta);

// Poisson link inverse: pmax(exp(eta), .Machine$double.eps)
NumericVector linkinv_poisson(const VectorXd &eta);

// Gaussian mu.eta: returns a vector of ones,
// analogous to rep.int(1, length(eta))
NumericVector mu_eta_gaussian(const VectorXd &eta);

// Binomial mu.eta: delegate to the exported logit_mu_eta_cpp function.
// It is assumed that logit_mu_eta_cpp is declared elsewhere and available.
NumericVector mu_eta_binomial(const VectorXd &eta);

// Poisson mu.eta: computes pmax(exp(eta), .Machine$double.eps)
NumericVector mu_eta_poisson(const VectorXd &eta);

List fastglm(const mat &x, const colvec &y, const string &family,
             Nullable<colvec> start = R_NilValue,
             Nullable<NumericVector> weights = R_NilValue,
             Nullable<NumericVector> offset = R_NilValue,
             Nullable<NumericVector> etastart = R_NilValue,
             Nullable<NumericVector> mustart = R_NilValue, int method = 0,
             double tol = 1e-8, int maxit = 100);

#endif  // FIT_GLM_H
