#include "fastcpd_classes.h"
#include "fastcpd_constants.h"

namespace fastcpd::classes {

List Fastcpd::negative_log_likelihood(
    mat data,
    Nullable<colvec> theta,
    double lambda,
    bool cv,
    Nullable<colvec> start
) {
  if (theta.isNull()) {
    return negative_log_likelihood_wo_theta(data, lambda, cv, start);
  } else {
    return List::create(
      Named("value") =
        negative_log_likelihood_wo_cv(data, as<colvec>(theta), lambda)
    );
  }
}

List Fastcpd::negative_log_likelihood_wo_theta(
    mat data,
    double lambda,
    bool cv,
    Nullable<colvec> start
) {
  if (family == "lasso" && cv) {
    Environment glmnet = Environment::namespace_env("glmnet"),
                stats = Environment::namespace_env("stats");
    Function cv_glmnet = glmnet["cv.glmnet"],
        predict_glmnet = glmnet["predict.glmnet"],
              deviance = stats["deviance"];
    List out = cv_glmnet(
      data.cols(1, data.n_cols - 1), data.col(0), Named("family") = "gaussian"
    );
    colvec index_vec = as<colvec>(out["index"]),
              values = as<colvec>(deviance(out["glmnet.fit"]));
    S4 out_coef = predict_glmnet(
      out["glmnet.fit"],
      Named("s") = out["lambda.1se"],
      Named("type") = "coefficients",
      Named("exact") = false
    );
    vec glmnet_i = as<vec>(out_coef.slot("i"));
    vec glmnet_x = as<vec>(out_coef.slot("x"));
    vec par = zeros(data.n_cols - 1);
    for (unsigned int i = 1; i < glmnet_i.n_elem; i++) {
      par(glmnet_i(i) - 1) = glmnet_x(i);
    }
    return List::create(
      Named("par") = par, Named("value") = values(index_vec(1) - 1)
    );
  } else if (family == "lasso" && !cv) {
    Environment stats = Environment::namespace_env("stats"),
              glmnet = Environment::namespace_env("glmnet");
    Function deviance = stats["deviance"], glmnet_ = glmnet["glmnet"],
      predict_glmnet = glmnet["predict.glmnet"];
    List out = glmnet_(
      data.cols(1, data.n_cols - 1), data.col(0),
      Named("family") = "gaussian", Named("lambda") = lambda
    );
    S4 out_par = out["beta"];
    vec par_i = as<vec>(out_par.slot("i"));
    vec par_x = as<vec>(out_par.slot("x"));
    vec par = zeros(data.n_cols - 1);
    for (unsigned int i = 0; i < par_i.n_elem; i++) {
      par(par_i(i)) = par_x(i);
    }
    double value = as<double>(deviance(out));
    vec fitted_values = as<vec>(
      predict_glmnet(out, data.cols(1, data.n_cols - 1), Named("s") = lambda)
    );
    vec residuals = data.col(0) - fitted_values;
    return List::create(Named("par") = par,
                  Named("value") = value / 2,
                  Named("residuals") = residuals);
  } else if (
    family == "binomial" || family == "poisson" || family == "gaussian"
  ) {
    vec y = data.col(0);
    Environment fastglm = Environment::namespace_env("fastglm");
    Function fastglm_ = fastglm["fastglm"];
    List out;
    if (start.isNull()) {
      mat x = data.cols(1, data.n_cols - 1);
      out = fastglm_(x, y, family);
    } else {
      colvec start_ = as<colvec>(start);
      mat x = data.cols(1, data.n_cols - 1);
      out = fastglm_(x, y, family, Named("start") = start_);
    }
    vec par = as<vec>(out["coefficients"]);
    vec residuals = as<vec>(out["residuals"]);
    double value = out["deviance"];
    return List::create(Named("par") = par,
                  Named("value") = value / 2,
                  Named("residuals") = residuals);
  } else if (family == "arma") {
    Environment stats = Environment::namespace_env("stats");
    Function arima = stats["arima"];
    List out = arima(
      Named("x") = data.col(0),
      Named("order") = Rcpp::NumericVector::create(order(0), 0, order(1)),
      Named("include.mean") = false
    );
    colvec par = zeros(sum(order) + 1);
    par.rows(0, sum(order) - 1) = as<colvec>(out["coef"]);
    par(sum(order)) = as<double>(out["sigma2"]);

    return List::create(Named("par") = par,
                  Named("value") = -as<double>(out["loglik"]),
                  Named("residuals") = as<vec>(out["residuals"]));
  } else if (family == "mean") {
    rowvec par = mean(data, 0);
    mat residuals = data.each_row() - par;
    return List::create(
      Named("par") = par,
      Named("value") = data.n_rows / 2.0 * (
        log(2.0 * M_PI) * data.n_cols + log_det_sympd(mean_data_cov) +
          arma::trace(
            solve(mean_data_cov, residuals.t() * residuals)
          ) / data.n_rows
      ),
      Named("residuals") = residuals
    );
  } else if (family == "variance") {
    mat residuals = data.each_row() - variance_data_mean;
    mat par = residuals.t() * residuals / (data.n_rows - 1);
    double value =
      data.n_cols * (log(2.0 * M_PI) + (data.n_rows - 1) / data.n_rows);
    if (data.n_rows >= data.n_cols) {
      value += log_det_sympd(par);
    }
    value *= data.n_rows / 2.0;
    return List::create(
      Named("par") = par,
      Named("value") = value,
      Named("residuals") = residuals
    );
  } else if (family == "meanvariance" || family == "mv") {
    mat covariance = arma::cov(data);

    double value =
      data.n_cols * (log(2.0 * M_PI) + (data.n_rows - 1) / data.n_rows);
    if (data.n_rows >= data.n_cols) {
      value += log_det_sympd(
        covariance + epsilon * eye<mat>(data.n_cols, data.n_cols)
      );
    }
    value *= data.n_rows / 2.0;

    colvec par = zeros(data.n_cols * data.n_cols + data.n_cols);
    par.rows(0, data.n_cols - 1) = mean(data, 0).t();
    par.rows(data.n_cols, par.n_rows - 1) =
      covariance.reshape(data.n_cols * data.n_cols, 1);
    mat residuals = data.each_row() - par.rows(0, data.n_cols - 1).t();

    return List::create(
      Named("par") = par,
      Named("value") = value,
      Named("residuals") = residuals
    );
  } else {
    // # nocov start
    stop("This branch should not be reached at functions.cc: 103.");
    // # nocov end
  }
}

double Fastcpd::negative_log_likelihood_wo_cv(
    mat data,
    colvec theta,
    double lambda
) {
  vec y = data.col(0);
  if (family == "lasso" || family == "gaussian") {
    // Calculate negative log likelihood in gaussian family
    double penalty = lambda * accu(arma::abs(theta));
    mat x = data.cols(1, data.n_cols - 1);
    return accu(arma::pow(y - x * theta, 2)) / 2 + penalty;
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
    colvec reversed_theta = arma::reverse(theta);
    if (data.n_rows < max(order) + 1) {
      return 0;
    }
    colvec variance_term = zeros(data.n_rows);
    for (unsigned int i = max(order); i < data.n_rows; i++) {
      variance_term(i) = data(i, 0) - arma::dot(
          reversed_theta.rows(order(1) + 1, sum(order)),
          data.rows(i - order(0), i - 1)
        ) - arma::dot(
          reversed_theta.rows(1, order(1)),
          variance_term.rows(i - order(1), i - 1)
        );
    }
    return (log(2.0 * M_PI) +
      log(theta(sum(order)))) * (data.n_rows - 2) / 2.0 +
      arma::dot(variance_term, variance_term) / 2.0 / theta(sum(order));
  } else {
    // # nocov start
    stop("This branch should not be reached at functions.cc: 153.");
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
    mat reversed_data = arma::reverse(data, 0);
    colvec reversed_theta = arma::reverse(theta);
    if (data.n_rows < max(order) + 1) {
      return arma::ones(theta.n_elem);
    }
    colvec variance_term = zeros(data.n_rows);
    for (unsigned int i = max(order); i < data.n_rows; i++) {
      variance_term(i) = data(i, 0) - arma::dot(
          reversed_theta.rows(order(1) + 1, sum(order)),
          data.rows(i - order(0), i - 1)
        ) - arma::dot(
          reversed_theta.rows(1, order(1)),
          variance_term.rows(i - order(1), i - 1)
        );
    }
    colvec reversed_variance_term = arma::reverse(variance_term);
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
    hessian = (x.t() * x) * std::min(as_scalar(prob), min_prob);
  } else if (family == "lasso" || family == "gaussian") {
    hessian = x.t() * x;
  } else if (family == "arma") {
    // TODO(doccstat): Maybe we can store all these computations
    mat reversed_data = arma::reverse(data, 0);
    colvec reversed_theta = arma::reverse(theta);
    if (data.n_rows < max(order) + 1) {
      return arma::eye(theta.n_elem, theta.n_elem);
    }
    colvec variance_term = zeros(data.n_rows);
    for (unsigned int i = max(order); i < data.n_rows; i++) {
      variance_term(i) = data(i, 0) - arma::dot(
          reversed_theta.rows(order(1) + 1, sum(order)),
          data.rows(i - order(0), i - 1)
        ) - arma::dot(
          reversed_theta.rows(1, order(1)),
          variance_term.rows(i - order(1), i - 1)
        );
    }
    colvec reversed_variance_term = arma::reverse(variance_term);
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
    mat reversed_coef_phi = arma::reverse(phi_coefficient, 0),
        reversed_coef_psi = arma::reverse(psi_coefficient, 0);
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
