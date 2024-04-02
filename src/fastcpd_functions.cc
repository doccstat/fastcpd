#include "fastcpd_functions.h"

using ::fastcpd::classes::CostResult;

namespace fastcpd::functions {

CostResult negative_log_likelihood_arma(
  const mat data,
  const colvec order
) {
  Environment stats = Environment::namespace_env("stats");
  Function arima = stats["arima"];
  List out = arima(
    Named("x") = data.col(0),
    Named("order") = NumericVector::create(order(0), 0, order(1)),
    Named("include.mean") = false
  );
  colvec par = zeros(sum(order) + 1);
  par.rows(0, sum(order) - 1) = as<colvec>(out["coef"]);
  par(sum(order)) = as<double>(out["sigma2"]);

  return {{par}, {as<vec>(out["residuals"])}, -as<double>(out["loglik"])};
}

CostResult negative_log_likelihood_glm(
  const mat data,
  Nullable<colvec> start,
  const std::string family
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
  return {{par}, {residuals}, value / 2};
}

CostResult negative_log_likelihood_lasso_cv(const mat data) {
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
  return {{par}, {mat()}, values(index_vec(1) - 1)};
}

CostResult negative_log_likelihood_lasso_wo_cv(
  const mat data,
  const double lambda
) {
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
  return {{par}, {residuals}, value / 2};
}

CostResult negative_log_likelihood_mean(
  const mat data,
  const mat variance_estimate
) {
  rowvec data_n = data.row(data.n_rows - 1);
  rowvec data_1 = data.row(0);
  unsigned int p = data.n_cols;
  return {
    {zeros<colvec>(p - 1)},  // # nocov
    {zeros<mat>(data.n_rows, p - 1)},
    std::log(2.0 * M_PI) * data.n_cols + log_det_sympd(variance_estimate) *
      data.n_rows / 2.0 + (
      data_n(p - 1) - data_1(p - 1) - dot(
        data_n.subvec(0, p - 2) - data_1.subvec(0, p - 2),
        data_n.subvec(0, p - 2) - data_1.subvec(0, p - 2)
      ) / (data.n_rows - 1)
    ) / 2.0
  };
}

CostResult negative_log_likelihood_meanvariance(
  const mat data,
  const double epsilon
) {
  mat covariance = cov(data);

  double value = data.n_rows * data.n_cols * (std::log(2.0 * M_PI) + 1) / 2.0;
  if (data.n_rows >= data.n_cols) {
    value += data.n_rows * log_det_sympd(
      covariance + epsilon * eye<mat>(data.n_cols, data.n_cols)
    ) / 2.0;
  }

  colvec par = zeros(data.n_cols * data.n_cols + data.n_cols);
  par.rows(0, data.n_cols - 1) = mean(data, 0).t();
  par.rows(data.n_cols, par.n_rows - 1) =
    covariance.reshape(data.n_cols * data.n_cols, 1);
  mat residuals = data.each_row() - par.rows(0, data.n_cols - 1).t();

  return {{par}, {residuals}, value};
}

CostResult negative_log_likelihood_mgaussian(
  const mat data,
  const unsigned int p_response,
  const mat variance_estimate
) {
  mat x = data.cols(p_response, data.n_cols - 1);
  mat y = data.cols(0, p_response - 1);
  mat x_t_x;

  if (data.n_rows <= data.n_cols - p_response + 1) {
    x_t_x = eye<mat>(data.n_cols - p_response, data.n_cols - p_response);
  } else {
    x_t_x = x.t() * x;
  }

  mat par = solve(x_t_x, x.t()) * y;
  DEBUG_RCOUT(par);
  mat residuals = y - x * par;
  double value =
    p_response * std::log(2.0 * M_PI) + log_det_sympd(variance_estimate);
  value *= data.n_rows;
  value += trace(solve(variance_estimate, residuals.t() * residuals));
  DEBUG_RCOUT(value);
  return {{par}, {residuals}, value / 2};
}

CostResult negative_log_likelihood_variance(const mat data) {
  rowvec data_n = data.row(data.n_rows - 1);
  rowvec data_1 = data.row(0);
  const unsigned int p = sqrt(data.n_cols);
  double det_value =
    det(arma::reshape(data_n - data_1, p, p) / (data.n_rows - 1));
  if (det_value <= 0) {
    det_value = 1e-10;
  }
  return {
    {zeros<mat>(p, p)},
    {zeros<mat>(data.n_rows, p)},
    (std::log(2.0 * M_PI) * p + p + log(det_value)) * (data.n_rows - 1) / 2.0
  };
}

}  // namespace fastcpd::functions