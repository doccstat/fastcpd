#ifndef FASTCPD_CLASSES_H_
#define FASTCPD_CLASSES_H_

#include <memory>

#include "fastcpd_types.h"
#include "RcppClock.h"
#include "RProgress.h"

namespace fastcpd::classes {

struct ColMat {
  mat data;

  operator colvec() const;
  operator mat() const;
  operator rowvec() const;
};

struct CostResult {
  ColMat par;
  ColMat residuals;
  double value;
};

struct CostFunction {
  const Function cost;
  const mat& data;

  CostFunction(const Function& cost, const mat& data);
  CostResult operator() (
      const unsigned int segment_start,
      const unsigned int segment_end,
      const Nullable<colvec>& theta,
      const double lambda,  // UNUSED
      const bool cv,  // UNUSED
      const Nullable<colvec>& start  // UNUSED
  ) const;
};

struct CostGradient {
  const Function cost_gradient;
  const mat& data;

  CostGradient(const Function& cost_gradient, const mat& data);
  const colvec operator() (
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  ) const;
};

struct CostHessian {
  const Function cost_hessian;
  const mat& data;

  CostHessian(const Function& cost_hessian, const mat& data);
  const mat operator() (
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  ) const;
};

class Fastcpd {
 public:
  Fastcpd(
    const double beta,
    Nullable<Function> cost,
    const string cost_adjustment,
    Nullable<Function> cost_gradient,
    Nullable<Function> cost_hessian,
    const bool cp_only,
    const unsigned int d,
    mat data,
    const double epsilon,
    const string family,
    Nullable<Function> k,
    colvec line_search,
    const colvec lower,
    const double momentum_coef,
    const colvec order,
    const int p,
    const unsigned int p_response,
    const double pruning_coef,
    const string r_clock,
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
  colvec cost_update_gradient(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  // Function to calculate the Hessian matrix at the current data.
  //
  // @param data A data frame containing the data to be segmented.
  // @param theta Estimated theta from the previous iteration.
  //
  // @return Hessian at the current data.
  mat cost_update_hessian(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  // Set \code{theta_sum} for a specific column.
  void create_theta_sum(const unsigned int col, colvec new_theta_sum);

  // Get the value of \code{theta_sum}.
  mat get_theta_sum();

  CostResult get_nll_wo_theta(
      const unsigned int segment_start,
      const unsigned int segment_end,
      double lambda,
      bool cv,
      Nullable<colvec> start
  );

  double get_nll_wo_cv(
      const unsigned int segment_start,
      const unsigned int segment_end,
      colvec theta,
      double lambda
  );

  List run();

  // Update \code{theta_sum} for a specific column by adding to that column.
  void update_theta_sum(const unsigned int col, colvec new_theta_sum);

 private:
  // `act_num` is used in Lasso and Gaussian families only.
  colvec act_num;

  // `beta` is the initial cost value.
  double beta;

  // `cost` is the cost function to be used.
  Nullable<Function> cost;

  // Adjustment to the cost function.
  const string cost_adjustment;

  // Cost function. If the cost function is provided in R, this will be a
  // wrapper of the R function.
  function<CostResult(
      const unsigned int segment_start,
      const unsigned int segment_end,
      Nullable<colvec> theta,
      const double lambda,
      const bool cv,
      Nullable<colvec> start
  )> cost_function_wrapper;

  // `cost_gradient` is the gradient of the cost function to be used.
  Nullable<Function> cost_gradient;

  // Gradient of the cost function. If the cost function is provided in R, this
  // will be a wrapper of the R function.
  function<colvec(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  )> cost_gradient_wrapper;

  // `cost_hessian` is the Hessian of the cost function to be used.
  Nullable<Function> cost_hessian;

  // Hessian of the cost function. If the cost function is provided in R, this
  // will be a wrapper of the R function.
  function<mat(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  )> cost_hessian_wrapper;

  const bool cp_only;

  // Dimension of the data.
  const unsigned int d;

  // `data` is the data set to be segmented.
  mat data;

  // The number of data points.
  const unsigned int data_n_rows;

  // The number of data columns.
  const unsigned int data_n_cols;

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

  const colvec order;

  // `p` is the number of parameters to be estimated.
  const unsigned int p;

  // Number of response variables in regression.
  const unsigned int p_response;

  const double pruning_coef;

  const string r_clock;
  const bool r_progress;

  Rcpp::Clock rClock;
  std::unique_ptr<RProgress::RProgress> rProgress;

  // `segment_count` is the number of segments for initial guess.
  const int segment_count;

  // `segment_indices` is the indices of the segments.
  colvec segment_indices;

  // Create a matrix to store the estimated coefficients in each segment,
  // where each row represents estimated coefficients for a segment.
  mat segment_theta_hat;

  // Matrix storing the warm starts.
  mat start;

  // `theta_hat` stores the estimated coefficients up to the current data point.
  mat theta_hat;

  // `theta_sum` stores the sum of estimated coefficients up to the current data
  // point.
  mat theta_sum;

  const double trim;

  // Upper bound of the parameters to be estimated during the optimization.
  const colvec upper;

  const double vanilla_percentage;

  const mat variance_estimate;

  const bool warm_start;

  mat zero_data;

  // Stop the clock and create an R object with `name`.
  void create_clock_in_r(const std::string name);

  void create_cost_function_wrapper(Nullable<Function> cost);
  void create_cost_gradient_wrapper(Nullable<Function> cost_gradient);
  void create_cost_hessian_wrapper(Nullable<Function> cost_hessian);

  // Initialize \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
  void create_gradients();

  // Initialize theta_hat_t_t to be the estimate in the segment.
  void create_segment_statistics();

  double get_cost_adjustment_value(const unsigned nrows);

  // Solve logistic/poisson regression using Gradient Descent Extension to the
  // multivariate case
  //
  // @param data A data frame containing the data to be segmented.
  // @param theta Estimate of the parameters. If null, the function will
  //   estimate the parameters.
  // @param lambda Lambda for L1 regularization. Only used for lasso.
  // @param cv Whether to perform cross-validation to find the best lambda.
  // @param start Starting point for the optimization for warm start.
  //   only used in mean change and lm.
  //
  // @return Negative log likelihood of the corresponding data with the given
  //   family.
  CostResult get_cost_result(
      const unsigned int segment_start,
      const unsigned int segment_end,
      Nullable<colvec> theta,
      const double lambda,
      const bool cv = false,
      Nullable<colvec> start = R_NilValue
  );

  List get_cp_set(const colvec raw_cp_set, const double lambda);

  double get_cval_for_r_t_set(
    const int tau,
    const unsigned int i,
    const int t,
    double lambda
  );

  double get_cval_pelt(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const unsigned int i,
    const double lambda
  );

  double get_cval_sen(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const unsigned int i,
    const double lambda
  );

  CostResult get_nll_arma(
    const unsigned int segment_start,
    const unsigned int segment_end
  );

  CostResult get_nll_glm(
    const unsigned int segment_start,
    const unsigned int segment_end,
    Nullable<colvec> start
  );

  CostResult get_nll_lasso_cv(
    const unsigned int segment_start,
    const unsigned int segment_end
  );

  CostResult get_nll_lasso_wo_cv(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const double lambda
  );

  CostResult get_nll_mean(
    const unsigned int segment_start,
    const unsigned int segment_end
  );

  CostResult get_nll_meanvariance(
    const unsigned int segment_start,
    const unsigned int segment_end
  );

  CostResult get_nll_mgaussian(
    const unsigned int segment_start,
    const unsigned int segment_end
  );

  CostResult get_nll_variance(
    const unsigned int segment_start,
    const unsigned int segment_end
  );

  // Update \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
  //
  // @param data_segment A data frame containing a segment of the data.
  //
  // @return A list containing new values of \code{theta_hat}, \code{theta_sum},
  //   and \code{hessian}.
  CostResult get_optimized_cost(
    const unsigned int segment_start,
    const unsigned int segment_end
  );

  void update_cost_parameters(
      const unsigned int t,
      const unsigned int tau,
      const unsigned int i,
      Function k,
      const double lambda,
      const colvec& line_search
  );

  void update_cost_parameters_step(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const int i,
    const int tau,
    const int j,
    const double lambda,
    const colvec& line_search
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
    const unsigned int segment_start,
    const unsigned int segment_end,
    const int tau,
    const int i,
    Function k,
    colvec momentum,
    const double lambda,
    const colvec& line_search
  );

  // Adjust cost value for MBIC and MDL.
  double update_cost_value(double value, const unsigned int nrows);

  colvec update_cp_set(const colvec raw_cp_set);

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

  // Start the clock tick for `name`.
  void update_r_clock_tick(const std::string name);

  // Stop the clock tick for `name`.
  void update_r_clock_tock(const std::string name);

  void update_r_progress_start();
  void update_r_progress_tick();

  void update_start(const unsigned int col, const colvec start_col);

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
};

}  // namespace fastcpd::classes

#endif  // FASTCPD_CLASSES_H_
