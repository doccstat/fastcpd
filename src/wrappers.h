#ifndef WRAPPERS_H
#define WRAPPERS_H

#include "RcppArmadillo.h"

using ::Rcpp::List;
using ::Rcpp::Nullable;

namespace fastcpd::wrappers {

class CostFunction {
 public:
  CostFunction(Rcpp::Function cost);

  List operator()(
      arma::mat data,
      Nullable<arma::colvec> theta,
      std::string family,  // UNUSED
      double lambda,  // UNUSED
      bool cv,  // UNUSED
      Nullable<arma::colvec> start  // UNUSED
  );

 private:
  Rcpp::Function cost;
};

class CostGradient {
 public:
  CostGradient(Rcpp::Function cost_gradient);

  arma::colvec operator()(
      arma::mat data,
      arma::colvec theta,
      std::string family  // UNUSED
  );

 private:
  Rcpp::Function cost_gradient;
};

class CostHessian {
 public:
  CostHessian(Rcpp::Function cost_hessian);

  arma::mat operator()(
      arma::mat data,
      arma::colvec theta,
      std::string family,  // UNUSED
      double min_prob  // UNUSED
  );

 private:
  Rcpp::Function cost_hessian;
};

}  // namespace fastcpd::wrappers

#endif  // WRAPPERS_H
