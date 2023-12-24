#ifndef FASTCPD_CLASSES_H_
#define FASTCPD_CLASSES_H_

#include "fastcpd_types.h"

namespace fastcpd::classes {

class Fastcpd {
 public:
  Fastcpd(
    mat data,
    const double beta,
    const int p,
    const colvec order,
    const string family,
    const double vanilla_percentage,
    const int segment_count,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double epsilon,
    const double min_prob,
    const double momentum_coef,
    const colvec lower,
    const colvec upper,
    const mat mean_data_cov,
    const unsigned int p_response
  );
  // Return `err_sd`.
  colvec get_err_sd();

  // Update `err_sd` for a specific segment.
  void update_err_sd(const unsigned int segment_index, const double err_var);

  // Segment the whole data set evenly based on the number of segments
  // specified in the `segment_count` parameter.
  void create_segment_indices();

  // Update \code{theta_hat} for a specific column.
  void update_theta_hat(const unsigned int col, colvec new_theta_hat);

  // Append a new column to \code{theta_hat}.
  void update_theta_hat(colvec new_theta_hat);

  // Prune the columns of \code{theta_hat}.
  void update_theta_hat(ucolvec pruned_left);

  // Get the value of \code{theta_sum}.
  mat get_theta_sum();

  // Set \code{theta_sum} for a specific column.
  void create_theta_sum(const unsigned int col, colvec new_theta_sum);

  // Update \code{theta_sum} for a specific column by adding to that column.
  void update_theta_sum(const unsigned int col, colvec new_theta_sum);

  // Append a new column to \code{theta_sum}.
  void update_theta_sum(colvec new_theta_sum);

  // Prune the columns of \code{theta_sum}.
  void update_theta_sum(ucolvec pruned_left);

  // Update \code{hessian} for a specific slice.
  void update_hessian(const unsigned int slice, mat new_hessian);

  // Append a new slice to \code{hessian}.
  void update_hessian(mat new_hessian);

  // Prune the slices of \code{hessian}.
  void update_hessian(ucolvec pruned_left);

  // Initialize theta_hat_t_t to be the estimate in the segment.
  void create_segment_statistics();

  double get_beta();

  // Update \code{momentum}.
  void update_momentum(colvec new_momentum);

  // Initialize \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
  void create_gradients();

  // Append new values to \code{fastcpd_parameters}.
  void update_fastcpd_parameters(const unsigned int t);

  void wrap_cost(Nullable<Function> cost);
  void wrap_cost_gradient(Nullable<Function> cost_gradient);
  void wrap_cost_hessian(Nullable<Function> cost_hessian);

  void cost_update(
      const unsigned int t,
      const unsigned int tau,
      const unsigned int i,
      Function k,
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
  // @param momentum_coef Momentum coefficient to be applied to the current
  //   momentum.
  // @param epsilon Epsilon to avoid numerical issues. Only used for binomial
  //   and poisson.
  // @param min_prob Minimum probability to avoid numerical issues.
  //   Only used for poisson.
  // @param winsorise_minval Minimum value to be winsorised. Only used for
  //   poisson.
  // @param winsorise_maxval Maximum value to be winsorised. Only used for
  //   poisson.
  // @param lambda Lambda for L1 regularization. Only used for lasso.
  // @param cost_gradient Gradient for custom cost function.
  // @param cost_hessian Hessian for custom cost function.
  // @param lower A vector containing the lower bounds for the parameters.
  // @param upper A vector containing the upper bounds for the parameters.
  // @param line_search A vector containing the line search coefficients.
  // @param order Order of the time series models.
  //
  // @return A list containing new values of \code{theta_hat}, \code{theta_sum},
  //   \code{hessian}, and \code{momentum}.
  List cost_update_steps(
    const mat data,
    const int tau,
    const int i,
    Function k,
    colvec momentum,
    const double lambda,
    colvec line_search
  );

  void cost_update_step(
    const mat data,
    const int i,
    const int tau,
    const int j,
    const double lambda,
    const colvec line_search
  );

  // `cost` is the cost function to be used.
  Nullable<Function> cost;

  // `cost_gradient` is the gradient of the cost function to be used.
  Nullable<Function> cost_gradient;

  // `cost_hessian` is the Hessian of the cost function to be used.
  Nullable<Function> cost_hessian;

  // Cost function. If the cost function is provided in R, this will be a
  // wrapper of the R function.
  function<List(
      mat data,
      Nullable<colvec> theta,
      double lambda,
      bool cv,
      Nullable<colvec> start
  )> cost_function_wrapper;

  // Gradient of the cost function. If the cost function is provided in R, this
  // will be a wrapper of the R function.
  function<colvec(mat data, colvec theta)> cost_gradient_wrapper;

  // Hessian of the cost function. If the cost function is provided in R, this
  // will be a wrapper of the R function.
  function<mat(mat data, colvec theta)> cost_hessian_wrapper;

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
  // @param mean_data_cov Covariance matrix of the data,
  //   only used in mean change.
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
  // @param family Family of the model.
  // @param min_prob Minimum probability to avoid numerical issues.
  // @param order Order of the time series models.
  //
  // @return Hessian at the current data.
  mat cost_update_hessian(mat data, colvec theta);

  // Update \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
  //
  // @param data_segment A data frame containing a segment of the data.
  // @param lambda Lambda for L1 regularization.
  //
  // @return A list containing new values of \code{theta_hat}, \code{theta_sum},
  //   and \code{hessian}.
  List cost_optim(
      const mat data_segment,
      const double lambda
  );

  Nullable<colvec> get_segment_theta_hat(const unsigned int t);

 private:
  // `data` is the data set to be segmented.
  mat data;

  // `beta` is the initial cost value.
  double beta;

  // `error_sd` is used in Gaussian family only.
  colvec err_sd;

  // `n` is the number of data points.
  int n;

  // `p` is the number of parameters to be estimated.
  const int p;

  const colvec order;

  // `family` is the family of the model.
  const string family;

  const double vanilla_percentage;

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

  const double min_prob;

  // `segment_indices` is the indices of the segments.
  colvec segment_indices;

  // `theta_hat` stores the estimated coefficients up to the current data point.
  mat theta_hat;

  // `theta_sum` stores the sum of estimated coefficients up to the current data
  // point.
  mat theta_sum;

  // `hessian` stores the Hessian matrix up to the current data point.
  cube hessian;

  // Create a matrix to store the estimated coefficients in each segment,
  // where each row represents estimated coefficients for a segment.
  mat segment_theta_hat;

  // `act_num` is used in Lasso and Gaussian families only.
  colvec act_num;

  // Momentum will be used in the update step if `momentum_coef` is not 0.
  colvec momentum;

  const double momentum_coef;

  const colvec lower;

  const colvec upper;

  const mat mean_data_cov;
  rowvec variance_data_mean;

  const unsigned int p_response;
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
