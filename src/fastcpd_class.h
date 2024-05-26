#ifndef FASTCPD_CLASS_H_
#define FASTCPD_CLASS_H_

#include <memory>

#include "fastcpd_test.h"
#include "RcppClock.h"
#include "RProgress.h"

using ::arma::cube;
using ::arma::ucolvec;
using ::Rcpp::Function;
using ::Rcpp::List;
using ::std::string;
using ::std::string_view;
using ::std::unique_ptr;

using ::fastcpd::test::FastcpdTest;

constexpr char kRProgress[] = "[:bar] :current/:total in :elapsed";

namespace fastcpd::classes {

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

  List run();

 private:
  // `act_num` is used in Lasso and Gaussian families only.
  colvec act_num;

  // `beta` is the initial cost value.
  double beta;

  // `cost` is the cost function to be used.
  unique_ptr<Function> cost;

  // Adjustment to the cost function.
  const string cost_adjustment;

  // `cost_gradient` is the gradient of the cost function to be used.
  unique_ptr<Function> cost_gradient;

  // `cost_hessian` is the Hessian of the cost function to be used.
  unique_ptr<Function> cost_hessian;

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

  // Function to calculate the gradient at the current data.
  //
  // @param data A data frame containing the data to be segmented.
  // @param theta Estimated theta from the previous iteration.
  // @param family Family of the model.
  // @param order Order of the time series models.
  //
  // @return Gradient at the current data.
  colvec (Fastcpd::*get_gradient)(
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
  mat (Fastcpd::*get_hessian)(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  CostResult (Fastcpd::*get_nll_pelt)(
      const unsigned int segment_start,
      const unsigned int segment_end,
      const double lambda,
      const bool cv,
      const Nullable<colvec>& start
  );

  double (Fastcpd::*get_nll_sen)(
      const unsigned int segment_start,
      const unsigned int segment_end,
      colvec theta,
      double lambda
  );

  // `hessian` stores the Hessian matrix up to the current data point.
  cube hessian;

  unique_ptr<Function> k;

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
  unique_ptr<RProgress::RProgress> rProgress;

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

  double** zero_data_c;

  // Stop the clock and create an R object with `name`.
  void create_clock_in_r(const std::string name);

  void create_gets(
    Nullable<Function>& cost,
    Nullable<Function>& cost_gradient,
    Nullable<Function>& cost_hessian
  );

  // Initialize \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
  void create_gradients();

  // Initialize theta_hat_t_t to be the estimate in the segment.
  void create_segment_statistics();

  // Set \code{theta_sum} for a specific column.
  void create_theta_sum(const unsigned int col, colvec new_theta_sum);

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

  colvec get_gradient_arma(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  colvec get_gradient_binomial(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  colvec get_gradient_custom(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  colvec get_gradient_lm(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  colvec get_gradient_ma(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  colvec get_gradient_poisson(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  mat get_hessian_arma(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  mat get_hessian_binomial(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  mat get_hessian_custom(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  mat get_hessian_lm(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  mat get_hessian_ma(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  mat get_hessian_poisson(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const colvec& theta
  );

  CostResult get_nll_pelt_arma(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const double lambda,
    const bool cv,
    const Nullable<colvec>& start
  );

  CostResult get_nll_pelt_custom(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const double lambda,
    const bool cv,
    const Nullable<colvec>& start
  );

  CostResult get_nll_pelt_glm(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const double lambda,
    const bool cv,
    const Nullable<colvec>& start
  );

  CostResult get_nll_pelt_lasso(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const double lambda,
    const bool cv,
    const Nullable<colvec>& start
  );
  CostResult get_nll_pelt_mean(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const double lambda,
    const bool cv,
    const Nullable<colvec>& start
  );

  CostResult get_nll_pelt_meanvariance(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const double lambda,
    const bool cv,
    const Nullable<colvec>& start
  );

  CostResult get_nll_pelt_mgaussian(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const double lambda,
    const bool cv,
    const Nullable<colvec>& start
  );

  CostResult get_nll_pelt_variance(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const double lambda,
    const bool cv,
    const Nullable<colvec>& start
  );

  double get_nll_sen_arma(
    const unsigned int segment_start,
    const unsigned int segment_end,
    colvec theta,
    double lambda
  );

  double get_nll_sen_binomial(
    const unsigned int segment_start,
    const unsigned int segment_end,
    colvec theta,
    double lambda
  );

  double get_nll_sen_custom(
    const unsigned int segment_start,
    const unsigned int segment_end,
    colvec theta,
    double lambda
  );

  double get_nll_sen_lm(
    const unsigned int segment_start,
    const unsigned int segment_end,
    colvec theta,
    double lambda
  );

  double get_nll_sen_ma(
    const unsigned int segment_start,
    const unsigned int segment_end,
    colvec theta,
    double lambda
  );

  double get_nll_sen_poisson(
    const unsigned int segment_start,
    const unsigned int segment_end,
    colvec theta,
    double lambda
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

  // Update \code{theta_sum} for a specific column by adding to that column.
  void update_theta_sum(const unsigned int col, colvec new_theta_sum);

  friend FastcpdTest;
};

}  // namespace fastcpd::classes

#endif  // FASTCPD_CLASS_H_
