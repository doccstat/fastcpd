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

bool valideta_gaussian(VectorXd const& eta);

bool valideta_binomial(VectorXd const& eta);

bool valideta_poisson(VectorXd const& eta);

bool validmu_gaussian(VectorXd const& mu);

bool validmu_binomial(VectorXd const& mu);

bool validmu_poisson(VectorXd const& mu);

// Gaussian deviance residuals: wt * ((y - mu)^2)
NumericVector dev_resids_gaussian(Map<VectorXd> const& y, VectorXd const& mu,
                                  Map<VectorXd> const& wt);

// Binomial deviance residuals: simply call the exported C function.
// Note: binomial_dev_resids_cpp is assumed to wrap the SEXP
// version defined elsewhere.
NumericVector dev_resids_binomial(Map<VectorXd> const& y, VectorXd const& mu,
                                  Map<VectorXd> const& wt);

// Poisson deviance residuals:
//   r = mu * wt,
//   for indices where y > 0, set r = wt * (y * log(y/mu) - (y - mu))
//   and return 2 * r.
NumericVector dev_resids_poisson(Map<VectorXd> const& y, VectorXd const& mu,
                                 Map<VectorXd> const& wt);

// Gaussian variance: always 1 for each element.
NumericVector var_gaussian(VectorXd const& mu);

// Binomial variance: mu * (1 - mu)
NumericVector var_binomial(VectorXd const& mu);

// Poisson variance: just mu
NumericVector var_poisson(VectorXd const& mu);

// Gaussian link inverse: identity (eta).
NumericVector linkinv_gaussian(VectorXd const& eta);

// Binomial link inverse: delegate to logit_linkinv_cpp.
NumericVector linkinv_binomial(VectorXd const& eta);

// Poisson link inverse: pmax(exp(eta), .Machine$double.eps)
NumericVector linkinv_poisson(VectorXd const& eta);

// Gaussian mu.eta: returns a vector of ones,
// analogous to rep.int(1, length(eta))
NumericVector mu_eta_gaussian(VectorXd const& eta);

// Binomial mu.eta: delegate to the exported logit_mu_eta_cpp function.
// It is assumed that logit_mu_eta_cpp is declared elsewhere and available.
NumericVector mu_eta_binomial(VectorXd const& eta);

// Poisson mu.eta: computes pmax(exp(eta), .Machine$double.eps)
NumericVector mu_eta_poisson(VectorXd const& eta);

List fastglm(mat const& x, colvec const& y, string const& family,
             std::optional<colvec> start = std::nullopt,
             Nullable<NumericVector> weights = R_NilValue,
             Nullable<NumericVector> offset = R_NilValue,
             Nullable<NumericVector> etastart = R_NilValue,
             Nullable<NumericVector> mustart = R_NilValue, int method = 0,
             double tol = 1e-8, int maxit = 100);

#endif  // FIT_GLM_H
