#include "functions.h"

namespace fastcpd::functions {

List negative_log_likelihood(
    mat data,
    Nullable<colvec> theta,
    string family,
    double lambda,
    bool cv,
    Nullable<colvec> start
) {
  if (theta.isNull()) {
    return negative_log_likelihood_wo_theta(data, family, lambda, cv, start);
  } else {
    return List::create(
      Named("value") = negative_log_likelihood_wo_cv(
        data, as<colvec>(theta), family, lambda, start
      )
    );
  }
}

List negative_log_likelihood_wo_theta(
    mat data,
    string family,
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
    mat x = data.cols(1, data.n_cols - 1);
    vec y = data.col(0);
    Environment fastglm = Environment::namespace_env("fastglm");
    Function fastglm_ = fastglm["fastglm"];
    List out;
    if (start.isNull()) {
      out = fastglm_(x, y, family);
    } else {
      colvec start_ = as<colvec>(start);
      out = fastglm_(x, y, family, Named("start") = start_);
    }
    vec par = as<vec>(out["coefficients"]);
    vec residuals = as<vec>(out["residuals"]);
    double value = out["deviance"];
    return List::create(Named("par") = par,
                  Named("value") = value / 2,
                  Named("residuals") = residuals);
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
    Nullable<colvec> start
) {
  vec y = data.col(0);
  mat x = data.cols(1, data.n_cols - 1);
  if (family == "lasso" || family == "gaussian") {
    // Calculate negative log likelihood in gaussian family
    double penalty = lambda * accu(arma::abs(theta));
    return accu(arma::pow(y - x * theta, 2)) / 2 + penalty;
  } else if (family == "binomial") {
    // Calculate negative log likelihood in binomial family
    colvec u = x * theta;
    return accu(-y % u + arma::log(1 + exp(u)));
  } else {
    // Calculate negative log likelihood in poisson family
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
  }
}

colvec cost_update_gradient(
    mat data,
    colvec theta,
    string family
) {
  rowvec new_data = data.row(data.n_rows - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double y = new_data(0);
  colvec gradient;
  if (family.compare("binomial") == 0) {
    gradient = - (y - 1 / (1 + exp(-as_scalar(x * theta)))) * x.t();
  } else if (family.compare("poisson") == 0) {
    gradient = - (y - exp(as_scalar(x * theta))) * x.t();
  } else {
    // `family` is either "lasso" or "gaussian".
    gradient = - (y - as_scalar(x * theta)) * x.t();
  }
  return gradient;
}

mat cost_update_hessian(
    mat data,
    colvec theta,
    string family,
    double min_prob
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
  } else {
    // `family` is either "lasso" or "gaussian".
    hessian = x.t() * x;
  }
  return hessian;
}

}
