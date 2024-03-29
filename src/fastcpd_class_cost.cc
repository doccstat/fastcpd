#include "fastcpd_classes.h"
#include "fastcpd_functions.h"

using ::fastcpd::functions::negative_log_likelihood_arma;
using ::fastcpd::functions::negative_log_likelihood_glm;
using ::fastcpd::functions::negative_log_likelihood_lasso_cv;
using ::fastcpd::functions::negative_log_likelihood_lasso_wo_cv;
using ::fastcpd::functions::negative_log_likelihood_mean;
using ::fastcpd::functions::negative_log_likelihood_meanvariance;
using ::fastcpd::functions::negative_log_likelihood_mgaussian;
using ::fastcpd::functions::negative_log_likelihood_variance;

namespace fastcpd::classes {

CostResult Fastcpd::negative_log_likelihood(
    mat data,
    Nullable<colvec> theta,
    double lambda,
    bool cv,
    Nullable<colvec> start
) {
  CostResult cost_result;
  if (theta.isNull()) {
    cost_result = negative_log_likelihood_wo_theta(data, lambda, cv, start);
  } else {
    cost_result = CostResult{
      colvec(),
      colvec(),
      negative_log_likelihood_wo_cv(data, as<colvec>(theta), lambda)
    };
  }
  cost_result.value = adjust_cost_value(cost_result.value, data.n_rows);
  return cost_result;
}

CostResult Fastcpd::negative_log_likelihood_wo_theta(
    mat data,
    double lambda,
    bool cv,
    Nullable<colvec> start
) {
  CostResult cost_result;
  if (family == "lasso" && cv) {
    cost_result = negative_log_likelihood_lasso_cv(data);
  } else if (family == "lasso" && !cv) {
    cost_result = negative_log_likelihood_lasso_wo_cv(data, lambda);
  } else if (
    family == "binomial" || family == "poisson" || family == "gaussian"
  ) {
    cost_result = negative_log_likelihood_glm(data, start, family);
  } else if (family == "arma") {
    cost_result = negative_log_likelihood_arma(data, order);
  } else if (family == "mean") {
    cost_result = negative_log_likelihood_mean(data, variance_estimate);
  } else if (family == "variance") {
    cost_result = negative_log_likelihood_variance(data, variance_data_mean);
  } else if (family == "meanvariance" || family == "mv") {
    cost_result = negative_log_likelihood_meanvariance(data, epsilon);
  } else if (family == "mgaussian") {
    cost_result = negative_log_likelihood_mgaussian(
      data, p_response, variance_estimate
    );
  } else {
    // # nocov start
    stop("This branch should not be reached at fastcpd_class_cost.cc: 193.");
    // # nocov end
  }
  return cost_result;
}

double Fastcpd::negative_log_likelihood_wo_cv(
    mat data,
    colvec theta,
    double lambda
) {
  vec y = data.col(0);
  if (family == "lasso" || family == "gaussian") {
    // Calculate negative log likelihood in gaussian family
    double penalty = lambda * accu(abs(theta));
    mat x = data.cols(1, data.n_cols - 1);
    return accu(square(y - x * theta)) / 2 + penalty;
  } else if (family == "binomial") {
    // Calculate negative log likelihood in binomial family
    mat x = data.cols(1, data.n_cols - 1);
    colvec u = x * theta;
    return accu(-y % u + arma::log(1 + exp(u)));
  } else if (family == "poisson") {
    mat x = data.cols(1, data.n_cols - 1);
    colvec u = x * theta;
    colvec y_factorial(y.n_elem);
    for (unsigned int i = 0; i < y.n_elem; i++) {
      double log_factorial = 0;
      for (int j = 1; j <= y(i); ++j) {
        log_factorial += std::log(j);
      }
      y_factorial(i) = log_factorial;
    }

    return accu(-y % u + exp(u) + y_factorial);
  } else if (family == "arma") {
    colvec reversed_theta = reverse(theta);
    if (data.n_rows < max(order) + 1) {
      return 0;
    }
    colvec variance_term = zeros(data.n_rows);
    for (unsigned int i = max(order); i < data.n_rows; i++) {
      variance_term(i) = data(i, 0) - dot(
          reversed_theta.rows(order(1) + 1, sum(order)),
          data.rows(i - order(0), i - 1)
        ) - dot(
          reversed_theta.rows(1, order(1)),
          variance_term.rows(i - order(1), i - 1)
        );
    }
    return (std::log(2.0 * M_PI) +
      std::log(theta(sum(order)))) * (data.n_rows - 2) / 2.0 +
      dot(variance_term, variance_term) / 2.0 / theta(sum(order));
  } else {
    // # nocov start
    stop("This branch should not be reached at functions.cc: 248.");
    // # nocov end
  }
}

colvec Fastcpd::cost_update_gradient(mat data, colvec theta) {
  rowvec new_data = data.row(data.n_rows - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double y = new_data(0);
  colvec gradient;
  if (family.compare("binomial") == 0) {
    gradient = - (y - 1 / (1 + exp(-as_scalar(x * theta)))) * x.t();
  } else if (family.compare("poisson") == 0) {
    gradient = - (y - exp(as_scalar(x * theta))) * x.t();
  } else if (family == "lasso" || family == "gaussian") {
    gradient = - (y - as_scalar(x * theta)) * x.t();
  } else if (family == "arma") {
    mat reversed_data = reverse(data, 0);
    colvec reversed_theta = reverse(theta);
    if (data.n_rows < max(order) + 1) {
      return ones(theta.n_elem);
    }
    colvec variance_term = zeros(data.n_rows);
    for (unsigned int i = max(order); i < data.n_rows; i++) {
      variance_term(i) = data(i, 0) - dot(
          reversed_theta.rows(order(1) + 1, sum(order)),
          data.rows(i - order(0), i - 1)
        ) - dot(
          reversed_theta.rows(1, order(1)),
          variance_term.rows(i - order(1), i - 1)
        );
    }
    colvec reversed_variance_term = reverse(variance_term);
    mat phi_coefficient = zeros(data.n_rows, order(0)),
        psi_coefficient = zeros(data.n_rows, order(1));
    for (unsigned int i = max(order); i < data.n_rows; i++) {
      phi_coefficient.row(i) = -reversed_data.rows(
        data.n_rows - i, data.n_rows - i + order(0) - 1
      ).t() - reversed_theta.rows(1, order(1)).t() *
      phi_coefficient.rows(i - order(1), i - 1);
    }
    for (unsigned int i = order(1); i < data.n_rows; i++) {
      psi_coefficient.row(i) = -reversed_variance_term.rows(
          data.n_rows - i, data.n_rows - i + order(1) - 1
        ).t() - reversed_theta.rows(1, order(1)).t() *
        psi_coefficient.rows(i - order(1), i - 1);
    }
    gradient = zeros(sum(order) + 1);
    gradient.rows(0, order(0) - 1) = phi_coefficient.row(data.n_rows - 1).t() *
      variance_term(data.n_rows - 1) / theta(sum(order));
    gradient.rows(order(0), sum(order) - 1) =
      psi_coefficient.row(data.n_rows - 1).t() *
      variance_term(data.n_rows - 1) / theta(sum(order));
    gradient(sum(order)) = 1.0 / 2.0 / theta(sum(order)) -
      std::pow(variance_term(data.n_rows - 1), 2) / 2.0 /
      std::pow(theta(sum(order)), 2);
  } else {
    // # nocov start
    stop("This branch should not be reached at functions.cc: 188.");
    // # nocov end
  }
  return gradient;
}

mat Fastcpd::cost_update_hessian(mat data, colvec theta) {
  rowvec new_data = data.row(data.n_rows - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  mat hessian;
  if (family.compare("binomial") == 0) {
    double prob = 1 / (1 + exp(-as_scalar(x * theta)));
    hessian = (x.t() * x) * as_scalar((1 - prob) * prob);
  } else if (family.compare("poisson") == 0) {
    double prob = exp(as_scalar(x * theta));
    // Prevent numerical issues if `prob` is too large.
    hessian = (x.t() * x) * std::min(as_scalar(prob), 1e10);
  } else if (family == "lasso" || family == "gaussian") {
    hessian = x.t() * x;
  } else if (family == "arma") {
    // TODO(doccstat): Maybe we can store all these computations
    mat reversed_data = reverse(data, 0);
    colvec reversed_theta = reverse(theta);
    if (data.n_rows < max(order) + 1) {
      return eye(theta.n_elem, theta.n_elem);
    }
    colvec variance_term = zeros(data.n_rows);
    for (unsigned int i = max(order); i < data.n_rows; i++) {
      variance_term(i) = data(i, 0) - dot(
          reversed_theta.rows(order(1) + 1, sum(order)),
          data.rows(i - order(0), i - 1)
        ) - dot(
          reversed_theta.rows(1, order(1)),
          variance_term.rows(i - order(1), i - 1)
        );
    }
    colvec reversed_variance_term = reverse(variance_term);
    mat phi_coefficient = zeros(data.n_rows, order(0)),
        psi_coefficient = zeros(data.n_rows, order(1));
    for (unsigned int i = max(order); i < data.n_rows; i++) {
      phi_coefficient.row(i) = -reversed_data.rows(
        data.n_rows - i, data.n_rows - i + order(0) - 1
      ).t() - reversed_theta.rows(1, order(1)).t() *
      phi_coefficient.rows(i - order(1), i - 1);
    }
    for (unsigned int i = order(1); i < data.n_rows; i++) {
      psi_coefficient.row(i) = -reversed_variance_term.rows(
          data.n_rows - i, data.n_rows - i + order(1) - 1
        ).t() - reversed_theta.rows(1, order(1)).t() *
        psi_coefficient.rows(i - order(1), i - 1);
    }
    mat reversed_coef_phi = reverse(phi_coefficient, 0),
        reversed_coef_psi = reverse(psi_coefficient, 0);
    cube phi_psi_coefficient = zeros(order(1), order(0), data.n_rows),
         psi_psi_coefficient = zeros(order(1), order(1), data.n_rows);
    for (unsigned int i = order(1); i < data.n_rows; i++) {
      mat phi_psi_coefficient_part = zeros(order(1), order(0)),
          psi_psi_coefficient_part = zeros(order(1), order(1));
      for (unsigned int j = 1; j <= order(1); j++) {
        phi_psi_coefficient_part +=
          phi_psi_coefficient.slice(i - j) * theta(order(0) - 1 + j);
      }
      phi_psi_coefficient.slice(i) = -reversed_coef_phi.rows(
        data.n_rows - i, data.n_rows - i + order(1) - 1
      ) - phi_psi_coefficient_part;
      for (unsigned int j = 1; j <= order(1); j++) {
        psi_psi_coefficient_part +=
          psi_psi_coefficient.slice(i - j) * theta(order(0) - 1 + j);
      }
      psi_psi_coefficient.slice(i) = -reversed_coef_psi.rows(
          data.n_rows - i, data.n_rows - i + order(1) - 1
        ) - reversed_coef_psi.rows(
          data.n_rows - i, data.n_rows - i + order(1) - 1
        ).t() - psi_psi_coefficient_part;
    }
    hessian = zeros(sum(order) + 1, sum(order) + 1);
    hessian.submat(0, 0, order(0) - 1, order(0) - 1) =
      phi_coefficient.row(data.n_rows - 1).t() *
      phi_coefficient.row(data.n_rows - 1) / theta(sum(order));
    hessian.submat(0, order(0), order(0) - 1, sum(order) - 1) = (
      phi_psi_coefficient.slice(data.n_rows - 1).t() *
        variance_term(data.n_rows - 1) +
        phi_coefficient.row(data.n_rows - 1).t() *
        psi_coefficient.row(data.n_rows - 1)
    ) / theta(sum(order));
    hessian.submat(order(0), 0, sum(order) - 1, order(0) - 1) =
      hessian.submat(0, order(0), order(0) - 1, sum(order) - 1).t();
    hessian.submat(0, sum(order), order(0) - 1, sum(order)) =
      -phi_coefficient.row(data.n_rows - 1).t() *
      variance_term(data.n_rows - 1) / theta(sum(order)) / theta(sum(order));
    hessian.submat(sum(order), 0, sum(order), order(0) - 1) =
      hessian.submat(0, sum(order), order(0) - 1, sum(order)).t();
    hessian.submat(order(0), order(0), sum(order) - 1, sum(order) - 1) = (
      psi_coefficient.row(data.n_rows - 1).t() *
      psi_coefficient.row(data.n_rows - 1) +
      psi_psi_coefficient.slice(data.n_rows - 1) *
      variance_term(data.n_rows - 1)
    ) / theta(sum(order));
    hessian.submat(order(0), sum(order), sum(order) - 1, sum(order)) =
      -psi_coefficient.row(data.n_rows - 1).t() *
      variance_term(data.n_rows - 1) / theta(sum(order)) / theta(sum(order));
    hessian.submat(sum(order), order(0), sum(order), sum(order) - 1) =
      hessian.submat(order(0), sum(order), sum(order) - 1, sum(order)).t();
    hessian(sum(order), sum(order)) =
      std::pow(variance_term(data.n_rows - 1), 2) /
      std::pow(theta(sum(order)), 3) -
      1.0 / 2.0 / std::pow(theta(sum(order)), 2);
  }
  return hessian;
}

}  // namespace fastcpd::classes
