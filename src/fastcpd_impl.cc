#include "fastcpd_impl.h"

#include "fastcpd.h"

// Implementation of the fastcpd algorithm.
//
// @param data A data frame containing the data to be segmented.
// @param beta Initial cost value.
// @param cost_adjustment Adjustment for the cost function.
// @param d Dimension of the data.
// @param segment_count Number of segments for initial guess.
// @param trim Trimming for the boundary change points.
// @param momentum_coef Momentum coefficient to be applied to each update.
// @param multiple_epochs_function Function on number of epochs in SGD.
// @param family Family of the models. Can be "binomial", "poisson", "lasso",
//   "lm" or "arma". If not provided, the user must specify the cost function
//   and its gradient (and Hessian).
// @param epsilon Epsilon to avoid numerical issues. Only used for binomial and
//   poisson.
// @param p Number of parameters to be estimated.
// @param order Order for time series models.
// @param cost Cost function to be used. If not specified, the default is
//   the negative log-likelihood for the corresponding family.
// @param cost_gradient Gradient for custom cost function.
// @param cost_hessian Hessian for custom cost function.
// @param cp_only Whether to return only the change points or with the cost
//   values for each segment. If family is not provided or set to be
//   "custom", this parameter will be set to be true.
// @param vanilla_percentage How many of the data should be processed through
//   vanilla PELT. Range should be between 0 and 1. If set to be 0, all data
//   will be processed through sequential gradient descnet. If set to be 1,
//   all data will be processed through vaniall PELT. If the cost function
//   have an explicit solution, i.e. does not depend on coefficients like
//   the mean change case, this parameter will be set to be 1.
// @param warm_start Whether to use warm start for the initial guess.
// @param lower A vector containing the lower bounds for the parameters.
// @param upper A vector containing the upper bounds for the parameters.
// @param line_search A vector containing the line search coefficients.
// @param variance_estimate Covariance matrix of the data, only used in mean
//   change and gaussian.
// @param p_response Dimension of the response, used with multivariate
//   response.
// @param pruning_coef The constant to satisfy the pruning condition.
// @param r_progress Whether to show progress bar.
//
// @return A list containing the change points and the cost values for each
//   segment.
// [[Rcpp::export]]
Rcpp::List fastcpd_impl(
    const arma::mat& data, const double beta,
    const std::string& cost_adjustment, const int segment_count,
    const double trim, const double momentum_coef,
    const Rcpp::Nullable<Rcpp::Function>& multiple_epochs_function,
    const std::string& family, const double epsilon, const int p,
    const arma::colvec& order, const Rcpp::Nullable<Rcpp::Function>& cost,
    const Rcpp::Nullable<Rcpp::Function>& cost_gradient,
    const Rcpp::Nullable<Rcpp::Function>& cost_hessian, const bool cp_only,
    const double vanilla_percentage, const bool warm_start,
    const arma::colvec& lower, const arma::colvec& upper,
    const arma::colvec& line_search, const arma::mat& variance_estimate,
    const unsigned int p_response, const double pruning_coef,
    const bool r_progress) {
  fastcpd::classes::Fastcpd fastcpd_class(
      beta, cost, cost_adjustment, cost_gradient, cost_hessian, cp_only, data,
      epsilon, family, multiple_epochs_function, line_search, lower,
      momentum_coef, order, p, p_response, pruning_coef, r_progress,
      segment_count, trim, upper, vanilla_percentage, variance_estimate,
      warm_start);
  return fastcpd_class.Run();
}
