#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "fastcpd_types.h"

namespace fastcpd::functions {

// Solve logistic/poisson regression using Gradient Descent Extension to the
// multivariate case
//
// @param data A data frame containing the data to be segmented.
// @param theta Estimate of the parameters. If null, the function will estimate
//   the parameters.
// @param family Family of the model.
// @param lambda Lambda for L1 regularization. Only used for lasso.
// @param cv Whether to perform cross-validation to find the best lambda.
// @param start Starting point for the optimization for warm start.
// @param order Order of the time series models.
// @param mean_data_cov Covariance matrix of the data, only used in mean change.
//
// @return Negative log likelihood of the corresponding data with the given
//   family.
List negative_log_likelihood(
    mat data,
    Nullable<colvec> theta,
    string family,
    double lambda,
    bool cv = false,
    Nullable<colvec> start = R_NilValue,
    const colvec order = colvec(1),
    const mat mean_data_cov = eye(1, 1)
);

List negative_log_likelihood_wo_theta(
    mat data,
    string family,
    double lambda,
    bool cv,
    Nullable<colvec> start,
    const colvec order,
    const mat mean_data_cov
);

double negative_log_likelihood_wo_cv(
    mat data,
    colvec theta,
    string family,
    double lambda,
    Nullable<colvec> start,
    const colvec order
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
    mat data,
    colvec theta,
    string family,
    const colvec order
);

// Function to calculate the Hessian matrix at the current data.
//
// @param data A data frame containing the data to be segmented.
// @param theta Estimated theta from the previous iteration.
// @param family Family of the model.
// @param min_prob Minimum probability to avoid numerical issues.
// @param order Order of the time series models.
//
// @return Hessian at the current data.
mat cost_update_hessian(
    mat data,
    colvec theta,
    string family,
    double min_prob,
    const colvec order
);

}  // namespace fastcpd::functions

#endif  // FUNCTIONS_H_
