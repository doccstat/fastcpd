#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <RcppArmadillo.h>

using ::Rcpp::Function;
using ::Rcpp::List;
using ::Rcpp::Nullable;

namespace fastcpd::parameters {

class FastcpdParameters {
 public:
  FastcpdParameters(
    arma::mat data,
    const double beta,
    const int p,
    const std::string family,
    const int segment_count,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double epsilon
  );
  // Return `err_sd`.
  arma::colvec get_err_sd();

  // Update `err_sd` for a specific segment.
  void update_err_sd(const unsigned int segment_index, const double err_var);

  // Segment the whole data set evenly based on the number of segments
  // specified in the `segment_count` parameter.
  void create_segment_indices();

  // Return the indices of the segments.
  arma::colvec get_segment_indices();

  // Get the value of \code{theta_hat}.
  arma::mat get_theta_hat();

  // Update \code{theta_hat} for a specific column.
  void update_theta_hat(const unsigned int col, arma::colvec new_theta_hat);

  // Append a new column to \code{theta_hat}.
  void update_theta_hat(arma::colvec new_theta_hat);

  // Prune the columns of \code{theta_hat}.
  void update_theta_hat(arma::ucolvec pruned_left);

  // Get the value of \code{theta_sum}.
  arma::mat get_theta_sum();

  // Set \code{theta_sum} for a specific column.
  void create_theta_sum(const unsigned int col, arma::colvec new_theta_sum);

  // Update \code{theta_sum} for a specific column by adding to that column.
  void update_theta_sum(const unsigned int col, arma::colvec new_theta_sum);

  // Append a new column to \code{theta_sum}.
  void update_theta_sum(arma::colvec new_theta_sum);

  // Prune the columns of \code{theta_sum}.
  void update_theta_sum(arma::ucolvec pruned_left);

  // Get the value of \code{hessian}.
  arma::cube get_hessian();

  // Update \code{hessian} for a specific slice.
  void update_hessian(const unsigned int slice, arma::mat new_hessian);

  // Append a new slice to \code{hessian}.
  void update_hessian(arma::mat new_hessian);

  // Prune the slices of \code{hessian}.
  void update_hessian(arma::ucolvec pruned_left);

  // Initialize theta_hat_t_t to be the estimate in the segment.
  void create_segment_statistics();

  // Adjust `beta` for Lasso and Gaussian families. This seems to be working
  // but there might be better choices.
  void update_beta();

  // Get the value of \code{momentum}.
  arma::colvec get_momentum();

  // Update \code{momentum}.
  void update_momentum(arma::colvec new_momentum);

  // Get the value of \code{segment_theta_hat}.
  arma::mat get_segment_theta_hat();

  // Get the value of \code{act_num}.
  arma::colvec get_act_num();

  // Initialize \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
  void create_gradients();

  // Append new values to \code{fastcpd_parameters}.
  void update_fastcpd_parameters(const unsigned int t);

  void wrap_cost(Nullable<Function> cost);
  void wrap_cost_gradient(Nullable<Function> cost_gradient);
  void wrap_cost_hessian(Nullable<Function> cost_hessian);
//   List cost_optim(
//       const arma::mat data_segment,
//       const double lambda,
//       const bool cv
//   );

  // `cost` is the cost function to be used.
  Nullable<Function> cost;

  // `cost_gradient` is the gradient of the cost function to be used.
  Nullable<Function> cost_gradient;

  // `cost_hessian` is the Hessian of the cost function to be used.
  Nullable<Function> cost_hessian;

  // Gradient of the cost function. If the cost function is provided in R, this
  // will be a wrapper of the R function.
  std::function<arma::colvec(
      arma::mat data,
      arma::colvec theta,
      std::string family
  )> cost_gradient_wrapper;

  // Hessian of the cost function. If the cost function is provided in R, this
  // will be a wrapper of the R function.
  std::function<arma::mat(
      arma::mat data,
      arma::colvec theta,
      std::string family,
      double min_prob
  )> cost_hessian_wrapper;

 private:
  // `data` is the data set to be segmented.
  arma::mat data;

  // `beta` is the initial cost value.
  double beta;

  // `error_sd` is used in Gaussian family only.
  arma::colvec err_sd;

  // `n` is the number of data points.
  int n;

  // `p` is the number of parameters to be estimated.
  const int p;

  // `family` is the family of the model.
  const std::string family;

  // `segment_count` is the number of segments for initial guess.
  const int segment_count;

  // `winsorise_minval` is the minimum value to be winsorised. Only used for
  // poisson.
  const double winsorise_minval;

  // `winsorise_maxval` is the maximum value to be winsorised. Only used for
  // poisson.
  const double winsorise_maxval;

  // `epsilon` is the epsilon to avoid numerical issues. Only used for binomial
  // and poisson.
  const double epsilon;

  // `segment_indices` is the indices of the segments.
  arma::colvec segment_indices;

  // `theta_hat` stores the estimated coefficients up to the current data point.
  arma::mat theta_hat;

  // `theta_sum` stores the sum of estimated coefficients up to the current data
  // point.
  arma::mat theta_sum;

  // `hessian` stores the Hessian matrix up to the current data point.
  arma::cube hessian;

  // Create a matrix to store the estimated coefficients in each segment,
  // where each row represents estimated coefficients for a segment.
  arma::mat segment_theta_hat;

  // `act_num` is used in Lasso and Gaussian families only.
  arma::colvec act_num;

  // Momentum will be used in the update step if `momentum_coef` is not 0.
  arma::colvec momentum;

  // Cost function. If the cost function is provided in R, this will be a
  // wrapper of the R function.
  std::function<List(
      arma::mat data,
      Nullable<arma::colvec> theta,
      std::string family,
      double lambda,
      bool cv,
      Nullable<arma::colvec> start
  )> cost_function_wrapper;

//   Function* winsorize;
//   void create_environment_functions();
};

}  // namespace fastcpd::parameters

#endif  // PARAMETERS_H
