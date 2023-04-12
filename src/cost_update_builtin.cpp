#include "fastcpd.h"

//' Update the cost values for the segmentation.
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
//'
//' @return A list containing new values of \code{theta_hat}, \code{theta_sum},
//'   \code{hessian}, and \code{momentum}.
// [[Rcpp::export]]
Rcpp::List cost_update_builtin(
    const arma::mat data,
    arma::mat theta_hat,
    arma::mat theta_sum,
    arma::cube hessian,
    const int tau,
    const int i,
    Rcpp::Function k,
    const std::string family,
    arma::colvec momentum,
    const double momentum_coef,
    const double epsilon,
    const double min_prob,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double lambda
) {
    // Get the hessian
    arma::mat hessian_i = hessian.slice(i - 1);
    hessian_i += cost_update_hessian(data, theta_hat.col(i - 1), family, min_prob);

    // Get the gradient
    arma::vec gradient = cost_update_gradient(data, theta_hat.col(i - 1), family);

    // Add epsilon to the diagonal for PSD hessian
    arma::mat hessian_psd = hessian_i + epsilon * arma::eye<arma::mat>(theta_hat.n_rows, theta_hat.n_rows);

    // Calculate momentum step
    arma::vec momentum_step = arma::solve(hessian_psd, gradient);
    momentum = momentum_coef * momentum - momentum_step;

    // Update theta_hat with momentum
    theta_hat.col(i - 1) += momentum;

    // Winsorize if family is Poisson
    if (family == "poisson") {
        Rcpp::Environment desc_tools("package:DescTools");
        Rcpp::Function winsorize = desc_tools["Winsorize"];
        Rcpp::NumericVector winsorize_result = winsorize(
            Rcpp::_["x"] = theta_hat.col(i - 1),
            Rcpp::_["minval"] = winsorise_minval,
            Rcpp::_["maxval"] = winsorise_maxval
        );
        theta_hat.col(i - 1) = arma::vec(winsorize_result.begin(), winsorize_result.size(), false);
    } else if (family == "lasso" || family == "gaussian") {
        // Update theta_hat with L1 penalty
        double hessian_norm = arma::norm(hessian_i, "fro");
        arma::vec normd = arma::abs(theta_hat.col(i - 1)) - lambda / hessian_norm;
        theta_hat.col(i - 1) = arma::sign(theta_hat.col(i - 1)) % arma::max(normd, 0.0);
    }

    for (int kk = 1; kk <= Rcpp::as<int>(k(data.n_rows - tau)); kk++) {
        for (int j = tau + 1; j <= data.n_rows; j++) {
            hessian_i += cost_update_hessian(
                data.rows(tau, j - 1), theta_hat.col(i - 1), family, min_prob
            );
            gradient = cost_update_gradient(
                data.rows(tau, j - 1), theta_hat.col(i - 1), family
            );
            hessian_psd = hessian_i + epsilon * arma::eye<arma::mat>(theta_hat.n_rows, theta_hat.n_rows);
            momentum_step = arma::solve(hessian_psd, gradient);
            momentum = momentum_coef * momentum - momentum_step;
            theta_hat.col(i - 1) += momentum;

            // Winsorize if family is Poisson
            if (family == "poisson") {
                Rcpp::Environment desc_tools("package:DescTools");
                Rcpp::Function winsorize = desc_tools["Winsorize"];
                Rcpp::NumericVector winsorize_result = winsorize(
                    Rcpp::_["x"] = theta_hat.col(i - 1),
                    Rcpp::_["minval"] = winsorise_minval,
                    Rcpp::_["maxval"] = winsorise_maxval
                );
                theta_hat.col(i - 1) = arma::vec(winsorize_result.begin(), winsorize_result.size(), false);
            } else if (family == "lasso" || family == "gaussian") {
                double hessian_norm = arma::norm(hessian_i, "fro");
                arma::vec normd = arma::abs(theta_hat.col(i - 1)) - lambda / hessian_norm;
                theta_hat.col(i - 1) = arma::sign(theta_hat.col(i - 1)) % arma::max(normd, 0.0);
            }
        }
    }

    theta_sum.col(i - 1) += theta_hat.col(i - 1);
    hessian.slice(i - 1) = hessian_i;
    return Rcpp::List::create(theta_hat.col(i - 1), theta_sum.col(i - 1), hessian.slice(i - 1), momentum);
}

