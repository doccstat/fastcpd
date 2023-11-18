#include "functions.h"

namespace fastcpd::functions {

List negative_log_likelihood(
    mat data,
    Nullable<colvec> theta,
    string family,
    double lambda,
    bool cv,
    Nullable<colvec> start,
    const colvec order
) {
  if (theta.isNull()) {
    return negative_log_likelihood_wo_theta(
      data, family, lambda, cv, start, order
    );
  } else {
    return List::create(
      Named("value") = negative_log_likelihood_wo_cv(
        data, as<colvec>(theta), family, lambda, start, order
      )
    );
  }
}

List negative_log_likelihood_wo_theta(
    mat data,
    string family,
    double lambda,
    bool cv,
    Nullable<colvec> start,
    const colvec order
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
    // TODO(doccstat): order and include.mean
    Environment stats = Environment::namespace_env("stats");
    Function arima = stats["arima"];
    List out = arima(
      Named("x") = data.col(0),
      Named("order") = Rcpp::IntegerVector::create(3, 0, 2),
      Named("include.mean") = false
    );
    colvec par = zeros(3 + 2 + 1);
    par.rows(0, 3 + 2 - 1) = as<colvec>(out["coef"]);
    par(3 + 2) = as<double>(out["sigma2"]);

    return List::create(Named("par") = par,
                  Named("value") = -as<double>(out["loglik"]),
                  Named("residuals") = as<vec>(out["residuals"]));
  } else {
    // # nocov start
    stop("This branch should not be reached at functions.cc: 103.");
    // # nocov end
  }
}

double negative_log_likelihood_wo_cv(
    mat data,
    colvec theta,
    string family,
    double lambda,
    Nullable<colvec> start,
    const colvec order
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
    // TODO(doccstat): order
    colvec reversed_theta = arma::reverse(theta);
    if (data.n_rows < std::max(3, 2) + 1) {
      return 0;
    }
    colvec variance_term = zeros(data.n_rows);
    for (unsigned int i = std::max(3, 2); i < data.n_rows; i++) {
      variance_term(i) = data(i, 0) -
        arma::dot(reversed_theta.rows(5 - 2, 5), data.rows(i - 3, i - 1)) -
        arma::dot(
          reversed_theta.rows(5 - 4, 5 - 3), variance_term.rows(i - 2, i - 1)
        );
    }
    return (log(2.0 * M_PI) + log(theta(5))) * (data.n_rows - 2) / 2.0 +
      arma::dot(variance_term, variance_term) / 2.0 / theta(5);
  } else {
    // # nocov start
    stop("This branch should not be reached at functions.cc: 153.");
    // # nocov end
  }
}

colvec cost_update_gradient(
    mat data,
    colvec theta,
    string family,
    const colvec order
) {
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
    // TODO(doccstat): order
    mat reversed_data = arma::reverse(data, 0);
    colvec reversed_theta = arma::reverse(theta);
    if (data.n_rows < std::max(3, 2) + 1) {
      return arma::ones(theta.n_elem);
    }
    colvec variance_term = zeros(data.n_rows);
    for (unsigned int i = std::max(3, 2); i < data.n_rows; i++) {
      variance_term(i) = data(i, 0) -
        arma::dot(reversed_theta.rows(5 - 2, 5), data.rows(i - 3, i - 1)) -
        arma::dot(
          reversed_theta.rows(5 - 4, 5 - 3), variance_term.rows(i - 2, i - 1)
        );
    }
    colvec reversed_variance_term = arma::reverse(variance_term);
    mat phi_coefficient = zeros(data.n_rows, 3),
        psi_coefficient = zeros(data.n_rows, 2);
    for (unsigned int i = std::max(3, 2); i < data.n_rows; i++) {
      phi_coefficient.row(i) =
        -reversed_data.rows(data.n_rows - i, data.n_rows - i + 2).t() -
        reversed_theta.rows(1, 2).t() * phi_coefficient.rows(i - 2, i - 1);
    }
    for (unsigned int i = 2; i < data.n_rows; i++) {
      psi_coefficient.row(i) =
        -reversed_variance_term.rows(data.n_rows - i, data.n_rows - i + 1).t() -
        reversed_theta.rows(1, 2).t() * psi_coefficient.rows(i - 2, i - 1);
    }
    gradient = zeros(3 + 2 + 1);
    gradient.rows(0, 2) = phi_coefficient.row(data.n_rows - 1).t() *
      variance_term(data.n_rows - 1) / theta(5);
    gradient.rows(3, 4) = psi_coefficient.row(data.n_rows - 1).t() *
      variance_term(data.n_rows - 1) / theta(5);
    gradient(5) = 1.0 / 2.0 / theta(5) -
      std::pow(variance_term(data.n_rows - 1), 2) / (2 * theta(5) * theta(5));
  } else {
    // # nocov start
    stop("This branch should not be reached at functions.cc: 188.");
    // # nocov end
  }
  return gradient;
}

mat cost_update_hessian(
    mat data,
    colvec theta,
    string family,
    double min_prob,
    const colvec order
) {
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
    // TODO(doccstat): order
    // TODO(doccstat): Maybe we can store all these computations
    mat reversed_data = arma::reverse(data, 0);
    colvec reversed_theta = arma::reverse(theta);
    if (data.n_rows < std::max(3, 2) + 1) {
      return arma::eye(theta.n_elem, theta.n_elem);
    }
    colvec variance_term = zeros(data.n_rows);
    for (unsigned int i = std::max(3, 2); i < data.n_rows; i++) {
      variance_term(i) = data(i, 0) -
        arma::dot(reversed_theta.rows(5 - 2, 5), data.rows(i - 3, i - 1)) -
        arma::dot(
          reversed_theta.rows(5 - 4, 5 - 3), variance_term.rows(i - 2, i - 1)
        );
    }
    colvec reversed_variance_term = arma::reverse(variance_term);
    mat phi_coefficient = zeros(data.n_rows, 3),
        psi_coefficient = zeros(data.n_rows, 2);
    for (unsigned int i = std::max(3, 2); i < data.n_rows; i++) {
      phi_coefficient.row(i) =
        -reversed_data.rows(data.n_rows - i, data.n_rows - i + 2).t() -
        reversed_theta.rows(1, 2).t() * phi_coefficient.rows(i - 2, i - 1);
    }
    for (unsigned int i = 2; i < data.n_rows; i++) {
      psi_coefficient.row(i) =
        -reversed_variance_term.rows(data.n_rows - i, data.n_rows - i + 1).t() -
        reversed_theta.rows(1, 2).t() * psi_coefficient.rows(i - 2, i - 1);
    }
    mat reversed_coef_phi = arma::reverse(phi_coefficient, 0),
        reversed_coef_psi = arma::reverse(psi_coefficient, 0);
    cube phi_psi_coefficient = zeros(2, 3, data.n_rows),
         psi_psi_coefficient = zeros(2, 2, data.n_rows);
    for (unsigned int i = 2; i < data.n_rows; i++) {
      mat phi_psi_coefficient_part = zeros(2, 3),
          psi_psi_coefficient_part = zeros(2, 2);
      for (unsigned int j = 1; j <= 2; j++) {
        phi_psi_coefficient_part +=
          phi_psi_coefficient.slice(i - j) * theta(2 + j);
      }
      phi_psi_coefficient.slice(i) =
        -reversed_coef_phi.rows(data.n_rows - i, data.n_rows - i + 1) -
          phi_psi_coefficient_part;
      for (unsigned int j = 1; j <= 2; j++) {
        psi_psi_coefficient_part +=
          psi_psi_coefficient.slice(i - j) * theta(2 + j);
      }
      psi_psi_coefficient.slice(i) =
        -reversed_coef_psi.rows(data.n_rows - i, data.n_rows - i + 1) -
        reversed_coef_psi.rows(data.n_rows - i, data.n_rows - i + 1).t() -
        psi_psi_coefficient_part;
    }
    hessian = zeros(3 + 2 + 1, 3 + 2 + 1);
    hessian.submat(0, 0, 2, 2) = phi_coefficient.row(data.n_rows - 1).t() *
      phi_coefficient.row(data.n_rows - 1) / theta(5);
    hessian.submat(0, 3, 2, 4) = (
      phi_psi_coefficient.slice(data.n_rows - 1).t() * variance_term(
        data.n_rows - 1
      ) + phi_coefficient.row(data.n_rows - 1).t() *
      psi_coefficient.row(data.n_rows - 1)
    ) / theta(5);
    hessian.submat(3, 0, 4, 2) = hessian.submat(0, 3, 2, 4).t();
    hessian.submat(0, 5, 2, 5) = -phi_coefficient.row(data.n_rows - 1).t() *
      variance_term(data.n_rows - 1) / theta(5) / theta(5);
    hessian.submat(5, 0, 5, 2) = hessian.submat(0, 5, 2, 5).t();
    hessian.submat(3, 3, 4, 4) = (
      psi_coefficient.row(data.n_rows - 1).t() *
      psi_coefficient.row(data.n_rows - 1) + psi_psi_coefficient.slice(
        data.n_rows - 1
      ) * variance_term(data.n_rows - 1)
    ) / theta(5);
    hessian.submat(3, 5, 4, 5) = -psi_coefficient.row(data.n_rows - 1).t() *
      variance_term(data.n_rows - 1) / theta(5) / theta(5);
    hessian.submat(5, 3, 5, 4) = hessian.submat(3, 5, 4, 5).t();
    hessian(5, 5) =
      std::pow(variance_term(data.n_rows - 1), 2) / std::pow(theta(5), 3) -
      1.0 / 2.0 / std::pow(theta(5), 2);
  }
  return hessian;
}

}
