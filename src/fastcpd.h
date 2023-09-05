#ifndef FASTCPD_H
#define FASTCPD_H

#include <RcppArmadillo.h>
#include <testthat.h>

using ::Rcpp::Function;
using ::Rcpp::List;
using ::Rcpp::Nullable;

// template <typename ReType, typename... Args>
// constexpr auto foo(ReType(*)(Args...)) noexcept {
//    return sizeof...(Args);
// }

class FastcpdParameters {
 public:
  FastcpdParameters(
    arma::mat data,
    const double beta,
    const int p,
    const std::string family,
    const int segment_count,
    Function cost,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double epsilon
  );
  // Create a matrix to store the estimated coefficients in each segment,
  // where each row represents estimated coefficients for a segment.
  arma::mat segment_theta_hat;

  // `error_sd` is used in Gaussian family only.
  arma::colvec err_sd;

  // `act_num` is used in Lasso and Gaussian families only.
  arma::colvec act_num;
  arma::mat theta_hat;
  arma::mat theta_sum;
  arma::cube hessian;

  // Momentum will be used in the update step if `momentum_coef` is not 0.
  arma::colvec momentum;
  void create_segment_indices();
  arma::colvec read_segment_indices();

  void update_theta_hat(const unsigned int col, arma::colvec new_theta_hat);
  void update_theta_hat(arma::colvec new_theta_hat);
  void create_theta_sum(const unsigned int col, arma::colvec new_theta_sum);
  void update_theta_sum(const unsigned int col, arma::colvec new_theta_sum);
  void update_theta_sum(arma::colvec new_theta_sum);
  void update_hessian(const unsigned int slice, arma::mat new_hessian);
  void update_hessian(arma::mat new_hessian);

  void update_theta_hat(arma::ucolvec pruned_left);
  void update_theta_sum(arma::ucolvec pruned_left);
  void update_hessian(arma::ucolvec pruned_left);

  // Initialize theta_hat_t_t to be the estimate in the segment.
  void create_segment_statistics();

  // Adjust `beta` for Lasso and Gaussian families. This seems to be working
  // but there might be better choices.
  void adjust_beta();

  // Initialize \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
  void create_gradients();

  // Append new values to \code{fastcpd_parameters}.
  void append_fastcpd_parameters(const unsigned int t);

 private:
  arma::mat data;
  double beta;
  int n;
  const int p;
  const std::string family;
  const int segment_count;
  Function cost;
  const double winsorise_minval;
  const double winsorise_maxval;
  const double epsilon;
  arma::colvec segment_indices;

//   Function* winsorize;
//   void create_environment_functions();
};

//' Solve logistic/poisson regression using Gradient Descent Extension to the
//' multivariate case
//' This function is not meant to be called directly by the user.
//'
//' @param data A data frame containing the data to be segmented.
//' @param theta Estimate of the parameters. If null, the function will estimate
//'   the parameters.
//' @param family Family of the model.
//' @param lambda Lambda for L1 regularization. Only used for lasso.
//' @param cv Whether to perform cross-validation to find the best lambda.
//' @param start Starting point for the optimization for warm start.
//' @keywords internal
//' @importFrom glmnet glmnet cv.glmnet predict.glmnet
//' @importFrom fastglm fastglm
//'
//' @noRd
//' @return Negative log likelihood of the corresponding data with the given
//'   family.
// [[Rcpp::export]]
List negative_log_likelihood(
    arma::mat data,
    Nullable<arma::colvec> theta,
    std::string family,
    double lambda,
    bool cv = false,
    Nullable<arma::colvec> start = R_NilValue
);

//' Function to calculate the gradient at the current data.
//' This function is not meant to be called directly by the user.
//'
//' @param data A data frame containing the data to be segmented.
//' @param theta Estimated theta from the previous iteration.
//' @param family Family of the model.
//' @keywords internal
//'
//' @noRd
//' @return Gradient at the current data.
// [[Rcpp::export]]
arma::colvec cost_update_gradient(
    arma::mat data,
    arma::colvec theta,
    std::string family
);

//' Function to calculate the Hessian matrix at the current data.
//' This function is not meant to be called directly by the user.
//'
//' @param data A data frame containing the data to be segmented.
//' @param theta Estimated theta from the previous iteration.
//' @param family Family of the model.
//' @param min_prob Minimum probability to avoid numerical issues.
//' @keywords internal
//'
//' @noRd
//' @return Hessian at the current data.
// [[Rcpp::export]]
arma::mat cost_update_hessian(
    arma::mat data,
    arma::colvec theta,
    std::string family,
    double min_prob
);

//' Update the cost values for the segmentation.
//' This function is not meant to be called directly by the user.
//'
//' @param data A data frame containing the data to be segmented.
//' @param theta_hat Estimated theta from the previous iteration.
//' @param theta_sum Sum of estimated theta from the previous iteration.
//' @param hessian Hessian matrix from the previous iteration.
//' @param tau Start of the current segment.
//' @param i Index of the current data in the whole data set.
//' @param k Number of epochs in SGD.
//' @param family Family of the model.
//' @param momentum Momentum from the previous iteration.
//' @param momentum_coef Momentum coefficient to be applied to the current
//'   momentum.
//' @param epsilon Epsilon to avoid numerical issues. Only used for binomial and
//'   poisson.
//' @param min_prob Minimum probability to avoid numerical issues. Only used for
//'   poisson.
//' @param winsorise_minval Minimum value to be winsorised. Only used for
//'   poisson.
//' @param winsorise_maxval Maximum value to be winsorised. Only used for
//'   poisson.
//' @param lambda Lambda for L1 regularization. Only used for lasso.
//' @param cost_gradient Gradient for custom cost function.
//' @param cost_hessian Hessian for custom cost function.
//' @keywords internal
//' @importFrom DescTools Winsorize
//'
//' @noRd
//' @return A list containing new values of \code{theta_hat}, \code{theta_sum},
//'   \code{hessian}, and \code{momentum}.
// [[Rcpp::export]]
List cost_update(
    const arma::mat data,
    arma::mat theta_hat,
    arma::mat theta_sum,
    arma::cube hessian,
    const int tau,
    const int i,
    Function k,
    const std::string family,
    arma::colvec momentum,
    const double momentum_coef,
    const double epsilon,
    const double min_prob,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double lambda,
    Function cost_gradient,
    Function cost_hessian
);

//' Update \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
//' This function is not meant to be called directly by the user.
//'
//' @param family Family of the model.
//' @param p Number of parameters.
//' @param data_segment A data frame containing a segment of the data.
//' @param cost Cost function.
//' @param lambda Lambda for L1 regularization.
//' @param cv Whether to perform cross-validation to find the best lambda.
//' @keywords internal
//'
//' @noRd
//' @return A list containing new values of \code{theta_hat}, \code{theta_sum},
//'   and \code{hessian}.
// [[Rcpp::export]]
List cost_optim(
    const std::string family,
    const int p,
    const arma::mat data_segment,
    Function cost,
    const double lambda,
    const bool cv
);

//' Implementation of the fastcpd algorithm.
//' This function is not meant to be called directly by the user.
//'
//' @param data A data frame containing the data to be segmented.
//' @param beta Initial cost value.
//' @param segment_count Number of segments for initial guess.
//' @param trim Trimming for the boundary change points.
//' @param momentum_coef Momentum coefficient to be applied to each update.
//' @param k Function on number of epochs in SGD.
//' @param family Family of the models. Can be "binomial", "poisson", "lasso" or
//'     "gaussian". If not provided, the user must specify the cost function and
//'     its gradient (and Hessian).
//' @param epsilon Epsilon to avoid numerical issues. Only used for binomial and
//'     poisson.
//' @param min_prob Minimum probability to avoid numerical issues. Only used for
//'     poisson.
//' @param winsorise_minval Minimum value to be winsorised. Only used for
//'     poisson.
//' @param winsorise_maxval Maximum value to be winsorised. Only used for
//'     poisson.
//' @param p Number of parameters to be estimated.
//' @param cost Cost function to be used. If not specified, the default is
//'     the negative log-likelihood for the corresponding family.
//' @param cost_gradient Gradient for custom cost function.
//' @param cost_hessian Hessian for custom cost function.
//' @param cp_only Whether to return only the change points or with the cost
//'     values for each segment. If family is not provided or set to be "custom",
//'     this parameter will be set to be true.
//' @param vanilla_percentage How many of the data should be processed through
//'     vanilla PELT. Range should be between 0 and 1. If set to be 0, all data
//'     will be processed through sequential gradient descnet. If set to be 1,
//'     all data will be processed through vaniall PELT. If the cost function
//'     have an explicit solution, i.e. does not depend on coefficients like
//'     the mean change case, this parameter will be set to be 1.
//' @keywords internal
//'
//' @noRd
//' @return A list containing the change points and the cost values for each
//'   segment.
// [[Rcpp::export]]
List fastcpd_impl(
    arma::mat data,
    double beta,
    const int segment_count,
    const double trim,
    const double momentum_coef,
    Function k,
    const std::string family,
    const double epsilon,
    const double min_prob,
    const double winsorise_minval,
    const double winsorise_maxval,
    const int p,
    Function cost,
    Function cost_gradient,
    Function cost_hessian,
    const bool cp_only,
    const double vanilla_percentage
);

#endif
