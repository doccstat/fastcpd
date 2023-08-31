#ifndef FASTCPD_H
#define FASTCPD_H

#include <RcppArmadillo.h>
#include <testthat.h>

using ::Rcpp::Function;
using ::Rcpp::List;
using ::Rcpp::Nullable;

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

//' Update the parameters related to fastcpd.
//' This function is not meant to be called directly by the user.
//'
//' @param fastcpd_parameters A list containing the parameters related to
//'   fastcpd.
//' @param data A data frame containing the data to be segmented.
//' @param t Current iteration.
//' @param i Index of the current data in the whole data set.
//' @param k Number of epochs in SGD.
//' @param tau Start of the current segment.
//' @param lambda Lambda for L1 regularization.
//' @param family Family of the model.
//' @param cost_gradient Gradient for custom cost function.
//' @param cost_hessian Hessian for custom cost function.
//' @param r_t_set Set of r_t values for the current iteration.
//' @param p Number of parameters.
//' @param momentum_coef Momentum coefficient to be applied to the current
//'   momentum.
//' @param min_prob Minimum probability to avoid numerical issues.
//' @param winsorise_minval Minimum value to be winsorised.
//' @param winsorise_maxval Maximum value to be winsorised.
//' @param epsilon Epsilon to avoid numerical issues.
//' @keywords internal
//'
//' @noRd
//' @return A list containing new values of \code{fastcpd_parameters}.
// [[Rcpp::export]]
List update_fastcpd_parameters(
    List fastcpd_parameters,
    arma::mat data,
    const int t,
    const int i,
    Function k,
    const int tau,
    const double lambda,
    const std::string family,
    Function cost_gradient,
    Function cost_hessian,
    arma::vec r_t_set,
    const int p,
    const double momentum_coef,
    const double min_prob,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double epsilon
);

//' Initialize \code{theta_hat}, \code{theta_sum}, and \code{hessian}.
//' This function is not meant to be called directly by the user.
//'
//' @param family Family of the model.
//' @param segment_theta_hat Estimated theta for the current segment.
//' @param data A data frame containing the data to be segmented.
//' @param p Number of parameters.
//' @param winsorise_minval Minimum value to be winsorised.
//' @param winsorise_maxval Maximum value to be winsorised.
//' @param epsilon Epsilon to avoid numerical issues.
//' @keywords internal
//'
//' @noRd
//' @return A list containing new values of \code{theta_hat}, \code{theta_sum},
//'   and \code{hessian}.
// [[Rcpp::export]]
List init_theta_hat_sum_hessian(
    const std::string family,
    const arma::mat segment_theta_hat,
    const arma::mat data,
    const int p,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double epsilon
);

#endif
