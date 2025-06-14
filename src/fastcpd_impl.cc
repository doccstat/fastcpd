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
    arma::mat const& data, double const beta,
    std::string const& cost_adjustment, int const segment_count,
    double const trim, double const momentum_coef,
    Rcpp::Nullable<Rcpp::Function> const& multiple_epochs_function,
    std::string const& family, double const epsilon, int const p,
    arma::colvec const& order, Rcpp::Nullable<Rcpp::Function> const& cost_pelt,
    Rcpp::Nullable<Rcpp::Function> const& cost_sen,
    Rcpp::Nullable<Rcpp::Function> const& cost_gradient,
    Rcpp::Nullable<Rcpp::Function> const& cost_hessian, bool const cp_only,
    double const vanilla_percentage, bool const warm_start,
    arma::colvec const& lower, arma::colvec const& upper,
    arma::colvec const& line_search, arma::mat const& variance_estimate,
    unsigned int const p_response, double const pruning_coef,
    bool const r_progress) {
  std::function<double(arma::mat)> cost_pelt_;
  if (family == "custom" && cost_pelt.isNotNull()) {
    // Capture the R function in a local Rcpp::Function object.
    Rcpp::Function rfun(cost_pelt);
    cost_pelt_ = [rfun](arma::mat data) -> double {
      // Call the R function and convert its result to double.
      return Rcpp::as<double>(rfun(data));
    };
  }
  std::function<double(arma::mat, arma::colvec)> cost_sen_;
  if (family == "custom" && cost_sen.isNotNull()) {
    // Capture the R function in a local Rcpp::Function object.
    Rcpp::Function rfun(cost_sen);
    cost_sen_ = [rfun](arma::mat data, arma::colvec theta) -> double {
      // Call the R function and convert its result to double.
      return Rcpp::as<double>(rfun(data, theta));
    };
  }
  std::optional<Rcpp::Function> cost = std::nullopt;
  if (family == "custom" && cost_pelt.isNotNull()) {
    Rcpp::Function rfun(cost_pelt);
    cost = rfun;
  } else if (family == "custom" && cost_sen.isNotNull()) {
    Rcpp::Function rfun(cost_sen);
    cost = rfun;
  }
  std::function<arma::colvec(arma::mat, arma::colvec)> cost_gradient_;
  if (family == "custom" && cost_gradient.isNotNull()) {
    // Capture the R function in a local Rcpp::Function object.
    Rcpp::Function rfun(cost_gradient);
    cost_gradient_ = [rfun](arma::mat data,
                            arma::colvec theta) -> arma::colvec {
      // Call the R function and convert its result to arma::colvec.
      return Rcpp::as<arma::colvec>(rfun(data, theta));
    };
  }
  std::function<arma::mat(arma::mat, arma::colvec)> cost_hessian_;
  if (family == "custom" && cost_hessian.isNotNull()) {
    // Capture the R function in a local Rcpp::Function object.
    Rcpp::Function rfun(cost_hessian);
    cost_hessian_ = [rfun](arma::mat data, arma::colvec theta) -> arma::mat {
      // Call the R function and convert its result to arma::mat.
      return Rcpp::as<arma::mat>(rfun(data, theta));
    };
  }
  std::function<unsigned int(unsigned int)> multiple_epochs_function_;
  if (multiple_epochs_function.isNotNull()) {
    // Capture the R function in a local Rcpp::Function object.
    Rcpp::Function rfun(multiple_epochs_function);
    multiple_epochs_function_ = [rfun](unsigned int i) -> unsigned int {
      // Call the R function and convert its result to unsigned int.
      return Rcpp::as<unsigned int>(rfun(i));
    };
  }
  fastcpd::classes::Fastcpd fastcpd_class(
      beta, cost, cost_pelt_, cost_sen_, cost_adjustment, cost_gradient_,
      cost_hessian_, cp_only, data, epsilon, family, multiple_epochs_function_,
      line_search, lower, momentum_coef, order, p, p_response, pruning_coef,
      r_progress, segment_count, trim, upper, vanilla_percentage,
      variance_estimate, warm_start);
  std::tuple<arma::colvec, arma::colvec, arma::colvec, arma::mat, arma::mat>
      result = fastcpd_class.Run();
  return Rcpp::List::create(Rcpp::Named("raw_cp_set") = std::get<0>(result),
                            Rcpp::Named("cp_set") = std::get<1>(result),
                            Rcpp::Named("cost_values") = std::get<2>(result),
                            Rcpp::Named("residual") = std::get<3>(result),
                            Rcpp::Named("thetas") = std::get<4>(result));
}
