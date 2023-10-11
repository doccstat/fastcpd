#ifndef FASTCPD_H
#define FASTCPD_H

#include <RcppArmadillo.h>

using ::Rcpp::Function;
using ::Rcpp::List;
using ::Rcpp::Nullable;

// Update the cost values for the segmentation.
//
// @param data A data frame containing the data to be segmented.
// @param theta_hat Estimated theta from the previous iteration.
// @param theta_sum Sum of estimated theta from the previous iteration.
// @param hessian Hessian matrix from the previous iteration.
// @param tau Start of the current segment.
// @param i Index of the current data in the whole data set.
// @param k Number of epochs in SGD.
// @param family Family of the model.
// @param momentum Momentum from the previous iteration.
// @param momentum_coef Momentum coefficient to be applied to the current
//   momentum.
// @param epsilon Epsilon to avoid numerical issues. Only used for binomial and
//   poisson.
// @param min_prob Minimum probability to avoid numerical issues. Only used for
//   poisson.
// @param winsorise_minval Minimum value to be winsorised. Only used for
//   poisson.
// @param winsorise_maxval Maximum value to be winsorised. Only used for
//   poisson.
// @param lambda Lambda for L1 regularization. Only used for lasso.
// @param cost_gradient Gradient for custom cost function.
// @param cost_hessian Hessian for custom cost function.
//
// @return A list containing new values of \code{theta_hat}, \code{theta_sum},
//   \code{hessian}, and \code{momentum}.
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
    std::function<arma::colvec(
        arma::mat data,
        arma::colvec theta,
        std::string family
    )> cost_gradient_wrapper,
    std::function<arma::mat(
        arma::mat data,
        arma::colvec theta,
        std::string family,
        double min_prob
    )> cost_hessian_wrapper
);

// Update \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
//
// @param family Family of the model.
// @param p Number of parameters.
// @param data_segment A data frame containing a segment of the data.
// @param cost Cost function.
// @param lambda Lambda for L1 regularization.
// @param cv Whether to perform cross-validation to find the best lambda.
//
// @return A list containing new values of \code{theta_hat}, \code{theta_sum},
//   and \code{hessian}.
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
//'   "gaussian". If not provided, the user must specify the cost function and
//'   its gradient (and Hessian).
//' @param epsilon Epsilon to avoid numerical issues. Only used for binomial and
//'   poisson.
//' @param min_prob Minimum probability to avoid numerical issues. Only used for
//'   poisson.
//' @param winsorise_minval Minimum value to be winsorised. Only used for
//'   poisson.
//' @param winsorise_maxval Maximum value to be winsorised. Only used for
//'   poisson.
//' @param p Number of parameters to be estimated.
//' @param cost Cost function to be used. If not specified, the default is
//'   the negative log-likelihood for the corresponding family.
//' @param cost_gradient Gradient for custom cost function.
//' @param cost_hessian Hessian for custom cost function.
//' @param cp_only Whether to return only the change points or with the cost
//'   values for each segment. If family is not provided or set to be
//'   "custom", this parameter will be set to be true.
//' @param vanilla_percentage How many of the data should be processed through
//'   vanilla PELT. Range should be between 0 and 1. If set to be 0, all data
//'   will be processed through sequential gradient descnet. If set to be 1,
//'   all data will be processed through vaniall PELT. If the cost function
//'   have an explicit solution, i.e. does not depend on coefficients like
//'   the mean change case, this parameter will be set to be 1.
//' @keywords internal
//' @importFrom DescTools Winsorize
//' @importFrom glmnet glmnet cv.glmnet predict.glmnet
//' @importFrom fastglm fastglm
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
    Nullable<Function> cost,
    Nullable<Function> cost_gradient,
    Nullable<Function> cost_hessian,
    const bool cp_only,
    const double vanilla_percentage
);

#endif  // FASTCPD_H
