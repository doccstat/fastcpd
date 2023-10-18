#ifndef COST_UPDATE_H_
#define COST_UPDATE_H_

#include "fastcpd_types.h"

namespace fastcpd::cost_update {

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
    const mat data,
    mat theta_hat,
    mat theta_sum,
    cube hessian,
    const int tau,
    const int i,
    Function k,
    const string family,
    colvec momentum,
    const double momentum_coef,
    const double epsilon,
    const double min_prob,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double lambda,
    function<colvec(mat data, colvec theta, string family)>
      cost_gradient_wrapper,
    function<mat(mat data, colvec theta, string family, double min_prob)>
      cost_hessian_wrapper
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
    const string family,
    const int p,
    const mat data_segment,
    Function cost,
    const double lambda,
    const bool cv
);

}  // namespace fastcpd::cost_update

#endif  // COST_UPDATE_H_
