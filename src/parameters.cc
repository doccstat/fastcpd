#include "constants.h"
#include "cost_update.h"
#include "fastcpd.h"
#include "functions.h"
#include "parameters.h"
#include "wrappers.h"

using ::fastcpd::cost_update::cost_optim;

namespace fastcpd::parameters {

FastcpdParameters::FastcpdParameters(
    mat data,
    const double beta,
    const int p,
    const string family,
    const int segment_count,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double epsilon
) : data(data),
    beta(beta),
    p(p),
    family(family),
    segment_count(segment_count),
    winsorise_minval(winsorise_minval),
    winsorise_maxval(winsorise_maxval),
    epsilon(epsilon) {
  n = data.n_rows;
  segment_indices = vec(n);
  segment_theta_hat = mat(segment_count, p);
  err_sd = vec(segment_count);
  act_num = vec(segment_count);
  theta_hat = mat(p, 1);
  theta_sum = mat(p, 1);
  hessian = cube(p, p, 1);
  momentum = vec(p);

  create_environment_functions();
}

void FastcpdParameters::create_environment_functions() {
  if (family == "poisson") {
    Environment desc_tools = Environment::namespace_env("DescTools");
    winsorize = as<Nullable<Function>>(desc_tools["Winsorize"]);
  } else {
    // TODO(doccstat): Store other environment functions.
  }
}

colvec FastcpdParameters::get_err_sd() {
  return err_sd;
}

void FastcpdParameters::update_err_sd(
  const unsigned int segment_index, const double err_var
) {
  err_sd(segment_index) = sqrt(err_var);
}

void FastcpdParameters::create_segment_indices() {
  const unsigned int segment_length = floor(n / segment_count);
  const int segment_remainder = n % segment_count;
  for (
    int segment_index = 1; segment_index <= segment_count; segment_index++
  ) {
    if (segment_index <= segment_remainder) {
      segment_indices(arma::span(
        (segment_index - 1) * (segment_length + 1),
        segment_index * (segment_length + 1)
      )).fill(segment_index);
    } else {
      segment_indices(arma::span(
        (segment_index - 1) * segment_length + segment_remainder,
        segment_index * segment_length + segment_remainder - 1
      )).fill(segment_index);
    }
  }
}

mat FastcpdParameters::get_theta_hat() {
  return theta_hat;
}

void FastcpdParameters::update_theta_hat(
  const unsigned int col, colvec new_theta_hat
) {
  theta_hat.col(col) = new_theta_hat;
}

void FastcpdParameters::update_theta_hat(colvec new_theta_hat) {
  theta_hat = arma::join_rows(theta_hat, new_theta_hat);
}

void FastcpdParameters::create_theta_sum(
  const unsigned int col, colvec new_theta_sum
) {
  theta_sum.col(col) = new_theta_sum;
}

mat FastcpdParameters::get_theta_sum() {
  return theta_sum;
}

void FastcpdParameters::update_theta_sum(
  const unsigned int col, colvec new_theta_sum
) {
  theta_sum.col(col) += new_theta_sum;
}

void FastcpdParameters::update_theta_sum(colvec new_theta_sum) {
  theta_sum = arma::join_rows(theta_sum, new_theta_sum);
}

cube FastcpdParameters::get_hessian() {
  return hessian;
}

void FastcpdParameters::update_hessian(
  const unsigned int slice, mat new_hessian
) {
  hessian.slice(slice) = new_hessian;
}

void FastcpdParameters::update_hessian(mat new_hessian) {
  hessian = arma::join_slices(hessian, new_hessian);
}

void FastcpdParameters::update_theta_hat(ucolvec pruned_left) {
  theta_hat = theta_hat.cols(pruned_left);
}

void FastcpdParameters::update_theta_sum(ucolvec pruned_left) {
  theta_sum = theta_sum.cols(pruned_left);
}

void FastcpdParameters::update_hessian(ucolvec pruned_left) {
  hessian = hessian.slices(pruned_left);
}

// TODO(doccstat): Use `segment_theta` as warm start.

void FastcpdParameters::create_segment_statistics() {
  for (
    int segment_index = 0; segment_index < segment_count; ++segment_index
  ) {
    ucolvec segment_indices_ =
        arma::find(segment_indices == segment_index + 1);
    mat data_segment = data.rows(segment_indices_);
    rowvec segment_theta;
    if (contain(CUSTOM_FAMILIES, family)) {
      segment_theta = as<rowvec>(
        cost_optim(family, p, data_segment, cost.get(), 0, false)["par"]
      );
    } else {
      segment_theta = as<rowvec>(
        cost_function_wrapper(
          data_segment, R_NilValue, family, 0, true, R_NilValue
        )["par"]
      );
    }

    // Initialize the estimated coefficients for each segment to be the
    // estimated coefficients in the segment.
    segment_theta_hat.row(segment_index) = segment_theta;
    if (family == "lasso" || family == "gaussian") {
      colvec segment_residual = data_segment.col(0) -
        data_segment.cols(1, data_segment.n_cols - 1) * segment_theta.t();
        double err_var = as_scalar(arma::mean(arma::pow(segment_residual, 2)));
        update_err_sd(segment_index, err_var);
        act_num(segment_index) = accu(arma::abs(segment_theta) > 0);
    }
  }
}

void FastcpdParameters::update_beta() {
  if (family == "lasso" || family == "gaussian") {
    beta = beta * (1 + mean(act_num));
  }
}

colvec FastcpdParameters::get_momentum() {
  return momentum;
}

void FastcpdParameters::update_momentum(colvec new_momentum) {
  momentum = new_momentum;
}

void FastcpdParameters::create_gradients() {
  if (family == "binomial") {
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    const double prob = 1 / (1 + exp(
      -as_scalar(theta_hat.t() * data.row(0).tail(data.n_cols - 1).t())
    ));
    hessian.slice(0) = (
      data.row(0).tail(data.n_cols - 1).t() * data.row(0).tail(data.n_cols - 1)
    ) * as_scalar(prob * (1 - prob));
  } else if (family == "poisson") {
    Function winsorize_non_null = winsorize.get();
    NumericVector winsorize_result = winsorize_non_null(
      Rcpp::_["x"] = segment_theta_hat.row(0).t(),
      Rcpp::_["minval"] = winsorise_minval,
      Rcpp::_["maxval"] = winsorise_maxval
    );
    theta_hat.col(0) = vec(
      winsorize_result.begin(), winsorize_result.size(), false
    );
    theta_sum.col(0) = vec(
      winsorize_result.begin(), winsorize_result.size(), false
    );
    hessian.slice(0) = (
      data.row(0).tail(data.n_cols - 1).t() * data.row(0).tail(data.n_cols - 1)
    ) * as_scalar(
      exp(theta_hat.t() * data.row(0).tail(data.n_cols - 1).t())
    );
  } else if (family == "lasso" || family == "gaussian") {
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    hessian.slice(0) = epsilon * eye<mat>(p, p) +
      data.row(0).tail(data.n_cols - 1).t() * data.row(0).tail(data.n_cols - 1);
  } else if (contain(CUSTOM_FAMILIES, family)) {
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    hessian.slice(0) = zeros<mat>(p, p);
  }
}

void FastcpdParameters::update_fastcpd_parameters(const unsigned int t) {
  // for tau = t-1
  rowvec new_data = data.row(t - 1).tail(data.n_cols - 1);
  const int segment_index = segment_indices(t - 1);
  rowvec cum_coef_add = segment_theta_hat.row(segment_index - 1),
             coef_add = segment_theta_hat.row(segment_index - 1);
  mat hessian_new;
  if (family == "binomial") {
    const double prob = 1 / (1 + exp(-as_scalar(coef_add * new_data.t())));
    hessian_new = (new_data.t() * new_data) * as_scalar(prob * (1 - prob));
  } else if (family == "poisson") {
    Function winsorize_non_null = winsorize.get();
    NumericVector winsorize_result = winsorize_non_null(
      Rcpp::_["x"] = coef_add,
      Rcpp::_["minval"] = winsorise_minval,
      Rcpp::_["maxval"] = winsorise_maxval
    );
    coef_add = as<rowvec>(winsorize_result);
    cum_coef_add = as<rowvec>(winsorize_result);
    hessian_new =
        (new_data.t() * new_data) * as_scalar(
          exp(coef_add * new_data.t())
        );
  } else if (family == "lasso" || family == "gaussian") {
    hessian_new = new_data.t() * new_data + epsilon * eye<mat>(p, p);
  } else if (contain(CUSTOM_FAMILIES, family)) {
    hessian_new = zeros<mat>(p, p);
  }

  update_theta_hat(coef_add.t());
  update_theta_sum(cum_coef_add.t());
  update_hessian(hessian_new);
}

void FastcpdParameters::wrap_cost(Nullable<Function> cost) {
  this->cost = cost;
  if (contain(FASTCPD_FAMILIES, family)) {
    cost_function_wrapper = &fastcpd::functions::negative_log_likelihood;
  } else if (cost.isNotNull()) {
    fastcpd::wrappers::CostFunction costFunction(cost.get());
    cost_function_wrapper = costFunction;
  } else if (cost.isNull()) {
    stop("cost function must be specified for custom family");
  } else {
    // This branch should not be reached.
  }
}

void FastcpdParameters::wrap_cost_gradient(Nullable<Function> cost_gradient) {
  this->cost_gradient = cost_gradient;
  if (contain(FASTCPD_FAMILIES, family)) {
    cost_gradient_wrapper = &fastcpd::functions::cost_update_gradient;
  } else if (cost_gradient.isNotNull()) {
    fastcpd::wrappers::CostGradient costGradient(cost_gradient.get());
    cost_gradient_wrapper = costGradient;
  } else if (cost_gradient.isNull()) {
    // `cost_gradient` can be `NULL` in the case of vanilla PELT.
  } else {
    // This branch should not be reached.
  }
}

void FastcpdParameters::wrap_cost_hessian(Nullable<Function> cost_hessian) {
  this->cost_hessian = cost_hessian;
  if (contain(FASTCPD_FAMILIES, family)) {
    cost_hessian_wrapper = &fastcpd::functions::cost_update_hessian;
  } else if (cost_hessian.isNotNull()) {
    fastcpd::wrappers::CostHessian costHessian(cost_hessian.get());
    cost_hessian_wrapper = costHessian;
  } else if (cost_hessian.isNull()) {
    // `cost_hessian` can be `NULL` in the case of vanilla PELT.
  } else {
    // This branch should not be reached.
  }
}

// List FastcpdParameters::cost_optim(
//     const mat data_segment,
//     const double lambda,
//     const bool cv
// ) {
//   if (family == "vanilla") {
//     return List::create(
//       Named("par") = R_NilValue,
//       Named("value") = cost_function_wrapper(
//         data_segment, R_NilValue, family, lambda, cv, R_NilValue
//       ),
//       Named("residuals") = R_NilValue
//     );
//   } else if (family == "custom" && p == 1) {
//     Environment stats = Environment::namespace_env("stats");
//     Function optim = stats["optim"];
//     List optim_result = optim(
//       Named("par") = 0,
//       Named("fn") = InternalFunction(
//         +[](double theta, mat data, Function cost) {
//           return cost(
//             Named("data") = data,
//             Named("theta") = log(theta / (1 - theta))
//           );
//         }
//       ),
//       Named("method") = "Brent",
//       Named("lower") = 0,
//       Named("upper") = 1,
//       Named("data") = data_segment,
//       Named("cost") = cost
//     );
//     return List::create(
//       Named("par") = log(
//         as<double>(optim_result["par"]) / (1 - as<double>(optim_result["par"]))
//       ),
//       Named("value") = exp(as<double>(optim_result["value"])) /
//         (1 + exp(as<double>(optim_result["value"]))),
//       Named("residuals") = R_NilValue
//     );
//   } else if (family == "custom" && p > 1) {
//     Environment stats = Environment::namespace_env("stats");
//     Function optim = stats["optim"];
//     List optim_result = optim(
//       Named("par") = zeros<vec>(p),
//       Named("fn") = cost,
//       Named("method") = "L-BFGS-B",
//       Named("data") = data_segment
//     );
//     return List::create(
//       Named("par") = optim_result["par"],
//       Named("value") = optim_result["value"],
//       Named("residuals") = R_NilValue
//     );
//   } else {
//     return cost(data_segment, R_NilValue, family, lambda, cv);
//   }
// }

}  // namespace fastcpd::parameters
