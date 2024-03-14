#ifndef FASTCPD_CLASSES_H_
#define FASTCPD_CLASSES_H_

#include "fastcpd_types.h"

namespace fastcpd::classes {

class Fastcpd {
 public:
  Fastcpd(
    const double beta,
    const double convexity_coef,
    Nullable<Function> cost,
    const string cost_adjustment,
    Nullable<Function> cost_gradient,
    Nullable<Function> cost_hessian,
    const bool cp_only,
    mat data,
    const double epsilon,
    const string family,
    Nullable<Function> k,
    colvec line_search,
    const colvec lower,
    const double momentum_coef,
    const colvec order,
    const int p,
    const bool pruning,
    const unsigned int p_response,
    const bool r_progress,
    const int segment_count,
    const double trim,
    const colvec upper,
    const double vanilla_percentage,
    const mat variance_estimate,
    const bool warm_start
  );

  // Function to calculate the gradient at the current data.
  //
  // @param data A data frame containing the data to be segmented.
  // @param theta Estimated theta from the previous iteration.
  // @param family Family of the model.
  // @param order Order of the time series models.
  //
  // @return Gradient at the current data.
  colvec cost_update_gradient(mat data, colvec theta);

  // Function to calculate the Hessian matrix at the current data.
  //
  // @param data A data frame containing the data to be segmented.
  // @param theta Estimated theta from the previous iteration.
  //
  // @return Hessian at the current data.
  mat cost_update_hessian(mat data, colvec theta);

  // Set \code{theta_sum} for a specific column.
  void create_theta_sum(const unsigned int col, colvec new_theta_sum);

  // Get the value of \code{theta_sum}.
  mat get_theta_sum();

  List negative_log_likelihood_wo_theta(
      mat data,
      double lambda,
      bool cv,
      Nullable<colvec> start
  );

  double negative_log_likelihood_wo_cv(
      mat data,
      colvec theta,
      double lambda
  );

  List run();

  // Update \code{theta_sum} for a specific column by adding to that column.
  void update_theta_sum(const unsigned int col, colvec new_theta_sum);

 private:

  // Adjust cost value for MBIC and MDL.
  double adjust_cost_value(double value, const unsigned int nrows);

  void create_cost_function_wrapper(Nullable<Function> cost);
  void create_cost_gradient_wrapper(Nullable<Function> cost_gradient);
  void create_cost_hessian_wrapper(Nullable<Function> cost_hessian);

  // Initialize \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
  void create_gradients();

  // Segment the whole data set evenly based on the number of segments
  // specified in the `segment_count` parameter.
  void create_segment_indices();

  // Initialize theta_hat_t_t to be the estimate in the segment.
  void create_segment_statistics();

  double get_cost_adjustment_value(const unsigned nrows);

  // Update \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
  //
  // @param data_segment A data frame containing a segment of the data.
  //
  // @return A list containing new values of \code{theta_hat}, \code{theta_sum},
  //   and \code{hessian}.
  List get_optimized_cost(const mat data_segment);

  // Solve logistic/poisson regression using Gradient Descent Extension to the
  // multivariate case
  //
  // @param data A data frame containing the data to be segmented.
  // @param theta Estimate of the parameters. If null, the function will
  //   estimate the parameters.
  // @param family Family of the model.
  // @param lambda Lambda for L1 regularization. Only used for lasso.
  // @param cv Whether to perform cross-validation to find the best lambda.
  // @param start Starting point for the optimization for warm start.
  // @param order Order of the time series models.
  // @param variance_estimate Covariance matrix of the data,
  //   only used in mean change and lm.
  //
  // @return Negative log likelihood of the corresponding data with the given
  //   family.
  List negative_log_likelihood(
      mat data,
      Nullable<colvec> theta,
      double lambda,
      bool cv = false,
      Nullable<colvec> start = R_NilValue
  );

  void update_cost_parameters(
      const unsigned int t,
      const unsigned int tau,
      const unsigned int i,
      Function k,
      const double lambda,
      const colvec line_search
  );

  void update_cost_parameters_step(
    const mat data,
    const int i,
    const int tau,
    const int j,
    const double lambda,
    const colvec line_search
  );

  // Update the cost values for the segmentation.
  //
  // @param data A data frame containing the data to be segmented.
  // @param tau Start of the current segment.
  // @param i Index of the current data in the whole data set.
  // @param k Number of epochs in SGD.
  // @param family Family of the model.
  // @param momentum Momentum from the previous iteration.
  // @param epsilon Epsilon to avoid numerical issues. Only used for binomial
  //   and poisson.
  // @param lambda Lambda for L1 regularization. Only used for lasso.
  // @param line_search A vector containing the line search coefficients.
  //
  // @return A list containing new values of \code{theta_hat}, \code{theta_sum},
  //   \code{hessian}, and \code{momentum}.
  List update_cost_parameters_steps(
    const mat data,
    const int tau,
    const int i,
    Function k,
    colvec momentum,
    const double lambda,
    colvec line_search
  );

  // Update `err_sd` for a specific segment.
  void update_err_sd(const unsigned int segment_index, const double err_var);

  // Append new values to \code{fastcpd_parameters}.
  void update_fastcpd_parameters(const unsigned int t);

  // Append a new slice to \code{hessian}.
  void update_hessian(mat new_hessian);

  // Prune the slices of \code{hessian}.
  void update_hessian(ucolvec pruned_left);

  // Update \code{hessian} for a specific slice.
  void update_hessian(const unsigned int slice, mat new_hessian);

  // Update \code{momentum}.
  void update_momentum(colvec new_momentum);

  // Append a new column to \code{theta_hat}.
  void update_theta_hat(colvec new_theta_hat);

  // Prune the columns of \code{theta_hat}.
  void update_theta_hat(ucolvec pruned_left);

  // Update \code{theta_hat} for a specific column.
  void update_theta_hat(const unsigned int col, colvec new_theta_hat);

  // Append a new column to \code{theta_sum}.
  void update_theta_sum(colvec new_theta_sum);

  // Prune the columns of \code{theta_sum}.
  void update_theta_sum(ucolvec pruned_left);

  // `act_num` is used in Lasso and Gaussian families only.
  colvec act_num;

  // `beta` is the initial cost value.
  double beta;

  const double convexity_coef;

  // `cost` is the cost function to be used.
  Nullable<Function> cost;

  // Adjustment to the cost function.
  const string cost_adjustment;

  // Cost function. If the cost function is provided in R, this will be a
  // wrapper of the R function.
  function<List(
      mat data,
      Nullable<colvec> theta,
      double lambda,
      bool cv,
      Nullable<colvec> start
  )> cost_function_wrapper;

  // `cost_gradient` is the gradient of the cost function to be used.
  Nullable<Function> cost_gradient;

  // Gradient of the cost function. If the cost function is provided in R, this
  // will be a wrapper of the R function.
  function<colvec(mat data, colvec theta)> cost_gradient_wrapper;

  // `cost_hessian` is the Hessian of the cost function to be used.
  Nullable<Function> cost_hessian;

  // Hessian of the cost function. If the cost function is provided in R, this
  // will be a wrapper of the R function.
  function<mat(mat data, colvec theta)> cost_hessian_wrapper;

  const bool cp_only;

  // `data` is the data set to be segmented.
  mat data;

  // `epsilon` is the epsilon to avoid numerical issues. Only used for binomial
  // and poisson.
  const double epsilon;

  // `error_sd` is used in Gaussian family only.
  colvec err_sd;

  // `family` is the family of the model.
  const string family;

  // `hessian` stores the Hessian matrix up to the current data point.
  cube hessian;

  Nullable<Function> k;

  colvec line_search;

  // Lower bound of the parameters to be estimated during the optimization.
  const colvec lower;

  // Momentum will be used in the update step if `momentum_coef` is not 0.
  colvec momentum;

  const double momentum_coef;

  // `n` is the number of data points.
  int n;

  const colvec order;

  // `p` is the number of parameters to be estimated.
  const int p;

  const bool pruning;

  // Number of response variables in regression.
  const unsigned int p_response;

  const bool r_progress;

  // `segment_count` is the number of segments for initial guess.
  const int segment_count;

  // `segment_indices` is the indices of the segments.
  colvec segment_indices;

  // Create a matrix to store the estimated coefficients in each segment,
  // where each row represents estimated coefficients for a segment.
  mat segment_theta_hat;

  // `theta_hat` stores the estimated coefficients up to the current data point.
  mat theta_hat;

  // `theta_sum` stores the sum of estimated coefficients up to the current data
  // point.
  mat theta_sum;

  const double trim;

  // Upper bound of the parameters to be estimated during the optimization.
  const colvec upper;

  const double vanilla_percentage;

  rowvec variance_data_mean;
  const mat variance_estimate;

  const bool warm_start;
};

class CostFunction {
 public:
  CostFunction(Function cost);

  List operator()(
      mat data,
      Nullable<colvec> theta,
      double lambda,  // UNUSED
      bool cv,  // UNUSED
      Nullable<colvec> start  // UNUSED
  );

 private:
  Function cost;
};

class CostGradient {
 public:
  CostGradient(Function cost_gradient);

  colvec operator()(mat data, colvec theta);

 private:
  Function cost_gradient;
};

class CostHessian {
 public:
  CostHessian(Function cost_hessian);

  mat operator()(mat data, colvec theta);

 private:
  Function cost_hessian;
};

}  // namespace fastcpd::parameters

#endif  // FASTCPD_CLASSES_H_
