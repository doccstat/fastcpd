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
//
// @return Negative log likelihood of the corresponding data with the given
//   family.
List negative_log_likelihood(
    mat data,
    Nullable<colvec> theta,
    string family,
    double lambda,
    bool cv = false,
    Nullable<colvec> start = R_NilValue
);

// Function to calculate the gradient at the current data.
//
// @param data A data frame containing the data to be segmented.
// @param theta Estimated theta from the previous iteration.
// @param family Family of the model.
//
// @return Gradient at the current data.
colvec cost_update_gradient(mat data, colvec theta, string family);

// Function to calculate the Hessian matrix at the current data.
//
// @param data A data frame containing the data to be segmented.
// @param theta Estimated theta from the previous iteration.
// @param family Family of the model.
// @param min_prob Minimum probability to avoid numerical issues.
//
// @return Hessian at the current data.
mat cost_update_hessian(
    mat data,
    colvec theta,
    string family,
    double min_prob
);

}  // namespace fastcpd::functions

#endif  // FUNCTIONS_H_
