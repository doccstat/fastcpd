#ifndef COST_UPDATE_H_
#define COST_UPDATE_H_

#include "fastcpd_types.h"

namespace fastcpd::cost_update {

// Update \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
//
// @param family Family of the model.
// @param p Number of parameters.
// @param data_segment A data frame containing a segment of the data.
// @param cost Cost function.
// @param lambda Lambda for L1 regularization.
// @param cv Whether to perform cross-validation to find the best lambda.
// @param lower A vector containing the lower bounds for the parameters.
// @param upper A vector containing the upper bounds for the parameters.
//
// @return A list containing new values of \code{theta_hat}, \code{theta_sum},
//   and \code{hessian}.
List cost_optim(
    const string family,
    const double vanilla_percentage,
    const int p,
    const mat data_segment,
    Function cost,
    const double lambda,
    const bool cv,
    const colvec lower,
    const colvec upper
);

}  // namespace fastcpd::cost_update

#endif  // COST_UPDATE_H_
