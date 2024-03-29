#include "fastcpd_functions.h"

using ::fastcpd::classes::CostResult;

namespace fastcpd::functions {

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
  return CostResult{par.t(), mat(), values(index_vec(1) - 1)};
}

CostResult negative_log_likelihood_mean(
  const mat data,
  const mat variance_estimate
) {
  rowvec par = mean(data, 0);
  return CostResult{
    par,
    data.each_row() - par,
    data.n_rows / 2.0 * (
      std::log(2.0 * M_PI) * data.n_cols + log_det_sympd(variance_estimate) +
        trace(solve(
          variance_estimate,
          (data.each_row() - par).t() * (data.each_row() - par)
        )) / data.n_rows
    )
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

  return CostResult{par.t(), residuals, value};
}

}  // namespace fastcpd::functions
