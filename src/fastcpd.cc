#include "fastcpd.h"

using ::Rcpp::as;
using ::Rcpp::Environment;
using ::Rcpp::Function;
using ::Rcpp::InternalFunction;
using ::Rcpp::List;
using ::Rcpp::Named;
using ::Rcpp::Nullable;
using ::Rcpp::NumericMatrix;
using ::Rcpp::NumericVector;
using ::Rcpp::S4;

// TODO(doccstat): `std::unordered_set` is not literal type and thus `constexpr`
// can not be used here. Blocked by `abseil` due to the non-header only issue.
const std::unordered_set<std::string> CUSTOM_FAMILIES = {
  "custom", "vanilla"
};

const std::unordered_set<std::string> FASTCPD_FAMILIES = {
  "gaussian", "binomial", "poisson", "lasso",
};

namespace fastcpd {

CostFunction::CostFunction(Rcpp::Function cost) : cost(cost) {}

List CostFunction::operator()(
    arma::mat data,
    Nullable<arma::colvec> theta,
    std::string family,  // UNUSED
    double lambda,  // UNUSED
    bool cv,  // UNUSED
    Nullable<arma::colvec> start  // UNUSED
) {
  return theta.isNull()? cost(data) : cost(data, theta);
}

CostGradient::CostGradient(Rcpp::Function cost_gradient) :
    cost_gradient(cost_gradient) {}

arma::colvec CostGradient::operator()(
    arma::mat data,
    arma::colvec theta,
    std::string family  // UNUSED
) {
  return as<arma::colvec>(cost_gradient(data, theta));
}

CostHessian::CostHessian(Rcpp::Function cost_hessian) :
    cost_hessian(cost_hessian) {}

arma::mat CostHessian::operator()(
    arma::mat data,
    arma::colvec theta,
    std::string family,  // UNUSED
    double min_prob  // UNUSED
) {
  return as<arma::mat>(cost_hessian(data, theta));
}

FastcpdParameters::FastcpdParameters(
    arma::mat data,
    const double beta,
    const int p,
    const std::string family,
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
  segment_indices = arma::vec(n);
  segment_theta_hat = arma::mat(segment_count, p);
  err_sd = arma::vec(segment_count);
  act_num = arma::vec(segment_count);
  theta_hat = arma::mat(p, 1);
  theta_sum = arma::mat(p, 1);
  hessian = arma::cube(p, p, 1);
  momentum = arma::vec(p);

  // create_environment_functions();
}

// void FastcpdParameters::create_environment_functions() {
//   if (family == "poisson") {
//     Environment desc_tools = Environment::namespace_env("DescTools");
//     *winsorize = desc_tools["Winsorize"];
//   }
// }

arma::colvec FastcpdParameters::get_err_sd() {
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

arma::vec FastcpdParameters::get_segment_indices() {
  return segment_indices;
}

arma::mat FastcpdParameters::get_theta_hat() {
  return theta_hat;
}

void FastcpdParameters::update_theta_hat(
  const unsigned int col, arma::colvec new_theta_hat
) {
  theta_hat.col(col) = new_theta_hat;
}

void FastcpdParameters::update_theta_hat(arma::colvec new_theta_hat) {
  theta_hat = arma::join_rows(theta_hat, new_theta_hat);
}

void FastcpdParameters::create_theta_sum(
  const unsigned int col, arma::colvec new_theta_sum
) {
  theta_sum.col(col) = new_theta_sum;
}

arma::mat FastcpdParameters::get_theta_sum() {
  return theta_sum;
}

void FastcpdParameters::update_theta_sum(
  const unsigned int col, arma::colvec new_theta_sum
) {
  theta_sum.col(col) += new_theta_sum;
}

void FastcpdParameters::update_theta_sum(arma::colvec new_theta_sum) {
  theta_sum = arma::join_rows(theta_sum, new_theta_sum);
}

arma::cube FastcpdParameters::get_hessian() {
  return hessian;
}

void FastcpdParameters::update_hessian(
  const unsigned int slice, arma::mat new_hessian
) {
  hessian.slice(slice) = new_hessian;
}

void FastcpdParameters::update_hessian(arma::mat new_hessian) {
  hessian = arma::join_slices(hessian, new_hessian);
}

void FastcpdParameters::update_theta_hat(arma::ucolvec pruned_left) {
  theta_hat = theta_hat.cols(pruned_left);
}

void FastcpdParameters::update_theta_sum(arma::ucolvec pruned_left) {
  theta_sum = theta_sum.cols(pruned_left);
}

void FastcpdParameters::update_hessian(arma::ucolvec pruned_left) {
  hessian = hessian.slices(pruned_left);
}

void FastcpdParameters::create_segment_statistics() {
  for (
    int segment_index = 0; segment_index < segment_count; ++segment_index
  ) {
    arma::ucolvec segment_indices_ =
        arma::find(segment_indices == segment_index + 1);
    arma::mat data_segment = data.rows(segment_indices_);
    arma::rowvec segment_theta;
    if (CUSTOM_FAMILIES.find(family) != CUSTOM_FAMILIES.end()) {
      segment_theta = as<arma::rowvec>(
        cost_optim(family, p, data_segment, cost.get(), 0, false)["par"]
      );
    } else {
      segment_theta = as<arma::rowvec>(negative_log_likelihood(
        data_segment, R_NilValue, family, 0, false
      )["par"]);
    }

    // Initialize the estimated coefficients for each segment to be the
    // estimated coefficients in the segment.
    segment_theta_hat.row(segment_index) = segment_theta;
    if (family == "lasso" || family == "gaussian") {
      arma::colvec segment_residual = data_segment.col(0) -
        data_segment.cols(1, data_segment.n_cols - 1) * segment_theta.t();
        double err_var =
            arma::as_scalar(arma::mean(arma::pow(segment_residual, 2)));
        update_err_sd(segment_index, err_var);
        act_num(segment_index) = arma::accu(arma::abs(segment_theta) > 0);
    }
  }
}

void FastcpdParameters::update_beta() {
  if (family == "lasso" || family == "gaussian") {
    beta = beta * (1 + mean(act_num));
  }
}

arma::colvec FastcpdParameters::get_momentum() {
  return momentum;
}

void FastcpdParameters::update_momentum(arma::colvec new_momentum) {
  momentum = new_momentum;
}

arma::mat FastcpdParameters::get_segment_theta_hat() {
  return segment_theta_hat;
}

arma::colvec FastcpdParameters::get_act_num() {
  return act_num;
}

void FastcpdParameters::create_gradients() {
  if (family == "binomial") {
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    const double prob = 1 / (1 + exp(
      -arma::as_scalar(theta_hat.t() * data.row(0).tail(data.n_cols - 1).t())
    ));
    hessian.slice(0) = (
      data.row(0).tail(data.n_cols - 1).t() * data.row(0).tail(data.n_cols - 1)
    ) * arma::as_scalar(prob * (1 - prob));
  } else if (family == "poisson") {
    Environment desc_tools = Environment::namespace_env("DescTools");
    Function winsorize = desc_tools["Winsorize"];
    NumericVector winsorize_result = winsorize(
      Rcpp::_["x"] = segment_theta_hat.row(0).t(),
      Rcpp::_["minval"] = winsorise_minval,
      Rcpp::_["maxval"] = winsorise_maxval
    );
    theta_hat.col(0) = arma::vec(
      winsorize_result.begin(), winsorize_result.size(), false
    );
    theta_sum.col(0) = arma::vec(
      winsorize_result.begin(), winsorize_result.size(), false
    );
    hessian.slice(0) = (
      data.row(0).tail(data.n_cols - 1).t() * data.row(0).tail(data.n_cols - 1)
    ) * arma::as_scalar(
      exp(theta_hat.t() * data.row(0).tail(data.n_cols - 1).t())
    );
  } else if (family == "lasso" || family == "gaussian") {
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    hessian.slice(0) = epsilon * arma::eye<arma::mat>(p, p) +
      data.row(0).tail(data.n_cols - 1).t() * data.row(0).tail(data.n_cols - 1);
  } else if (CUSTOM_FAMILIES.find(family) != CUSTOM_FAMILIES.end()) {
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    hessian.slice(0) = arma::zeros<arma::mat>(p, p);
  }
}

void FastcpdParameters::update_fastcpd_parameters(const unsigned int t) {
  // for tau = t-1
  arma::rowvec new_data = data.row(t - 1).tail(data.n_cols - 1);
  const int segment_index = segment_indices(t - 1);
  arma::rowvec cum_coef_add = segment_theta_hat.row(segment_index - 1),
                    coef_add = segment_theta_hat.row(segment_index - 1);
  arma::mat hessian_new;
  if (family == "binomial") {
    const double prob =
        1 / (1 + exp(-arma::as_scalar(coef_add * new_data.t())));
    hessian_new =
        (new_data.t() * new_data) * arma::as_scalar(prob * (1 - prob));
  } else if (family == "poisson") {
    Environment desc_tools = Environment::namespace_env("DescTools");
    Function winsorize = desc_tools["Winsorize"];
    NumericVector winsorize_result = winsorize(
      Rcpp::_["x"] = coef_add,
      Rcpp::_["minval"] = winsorise_minval,
      Rcpp::_["maxval"] = winsorise_maxval
    );
    coef_add = as<arma::rowvec>(winsorize_result);
    cum_coef_add = as<arma::rowvec>(winsorize_result);
    hessian_new =
        (new_data.t() * new_data) * arma::as_scalar(
          exp(coef_add * new_data.t())
        );
  } else if (family == "lasso" || family == "gaussian") {
    hessian_new =
        new_data.t() * new_data + epsilon * arma::eye<arma::mat>(p, p);
  } else if (CUSTOM_FAMILIES.find(family) != CUSTOM_FAMILIES.end()) {
    hessian_new = arma::zeros<arma::mat>(p, p);
  }

  update_theta_hat(coef_add.t());
  update_theta_sum(cum_coef_add.t());
  update_hessian(hessian_new);
}

void FastcpdParameters::wrap_cost(Nullable<Function> cost) {
  this->cost = cost;
  if (FASTCPD_FAMILIES.find(family) != FASTCPD_FAMILIES.end()) {
    cost_function_wrapper = &negative_log_likelihood;
  } else if (cost.isNotNull()) {
    CostFunction costFunction(cost.get());
    cost_function_wrapper = costFunction;
  } else if (cost.isNull()) {
    Rcpp::stop("cost function must be specified for custom family");
  } else {
    // NOTEST
    Rcpp::stop("Internal error, this branch should not be reached.");
  }
}

void FastcpdParameters::wrap_cost_gradient(Nullable<Function> cost_gradient) {
  this->cost_gradient = cost_gradient;
  if (FASTCPD_FAMILIES.find(family) != FASTCPD_FAMILIES.end()) {
    cost_gradient_wrapper = &cost_update_gradient;
  } else if (cost_gradient.isNotNull()) {
    CostGradient costGradient(cost_gradient.get());
    cost_gradient_wrapper = costGradient;
  } else if (cost_gradient.isNull()) {
    // `cost_gradient` can be `NULL` in the case of vanilla PELT.
  } else {
    // NOTEST
    Rcpp::stop("Internal error, this branch should not be reached.");
  }
}

void FastcpdParameters::wrap_cost_hessian(Nullable<Function> cost_hessian) {
  this->cost_hessian = cost_hessian;
  if (FASTCPD_FAMILIES.find(family) != FASTCPD_FAMILIES.end()) {
    cost_hessian_wrapper = &cost_update_hessian;
  } else if (cost_hessian.isNotNull()) {
    CostHessian costHessian(cost_hessian.get());
    cost_hessian_wrapper = costHessian;
  } else if (cost_hessian.isNull()) {
    // `cost_hessian` can be `NULL` in the case of vanilla PELT.
  } else {
    // NOTEST
    Rcpp::stop("Internal error, this branch should not be reached.");
  }
}

// List FastcpdParameters::cost_optim(
//     const arma::mat data_segment,
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
//         +[](double theta, arma::mat data, Function cost) {
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
//       Named("par") = arma::zeros<arma::vec>(p),
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

}  // namespace fastcpd

List negative_log_likelihood(
    arma::mat data,
    Nullable<arma::colvec> theta,
    std::string family,
    double lambda,
    bool cv,
    Nullable<arma::colvec> start
) {
  if (theta.isNull() && family == "lasso" && cv) {
    Environment glmnet = Environment::namespace_env("glmnet");
    Function cv_glmnet = glmnet["cv.glmnet"],
        predict_glmnet = glmnet["predict.glmnet"];
    List out = cv_glmnet(
      data.cols(1, data.n_cols - 1), data.col(0), Named("family") = "gaussian"
    );
    S4 out_coef = predict_glmnet(
      out["glmnet.fit"],
      Named("s") = out["lambda.1se"],
      Named("type") = "coefficients",
      Named("exact") = false
    );
    arma::vec glmnet_i = as<arma::vec>(out_coef.slot("i"));
    arma::vec glmnet_x = as<arma::vec>(out_coef.slot("x"));
    arma::vec par = arma::zeros(data.n_cols - 1);
    for (unsigned int i = 1; i < glmnet_i.n_elem; i++) {
      par(glmnet_i(i) - 1) = glmnet_x(i);
    }
    return List::create(Named("par") = par,
                  Named("value") = R_NilValue);
  } else if (theta.isNull() && family == "lasso" && !cv) {
    Environment stats = Environment::namespace_env("stats"),
               glmnet = Environment::namespace_env("glmnet");
    Function deviance = stats["deviance"], glmnet_ = glmnet["glmnet"],
       predict_glmnet = glmnet["predict.glmnet"];
    List out = glmnet_(
      data.cols(1, data.n_cols - 1), data.col(0),
      Named("family") = "gaussian", Named("lambda") = lambda
    );
    S4 out_par = out["beta"];
    arma::vec par_i = as<arma::vec>(out_par.slot("i"));
    arma::vec par_x = as<arma::vec>(out_par.slot("x"));
    arma::vec par = arma::zeros(data.n_cols - 1);
    for (unsigned int i = 0; i < par_i.n_elem; i++) {
      par(par_i(i)) = par_x(i);
    }
    double value = as<double>(deviance(out));
    arma::vec fitted_values = as<arma::vec>(
      predict_glmnet(out, data.cols(1, data.n_cols - 1), Named("s") = lambda)
    );
    arma::vec residuals = data.col(0) - fitted_values;
    return List::create(Named("par") = par,
                  Named("value") = value / 2,
                  Named("residuals") = residuals);
  } else if (theta.isNull()) {
    // Estimate theta in binomial/poisson/gaussian family
    arma::mat x = data.cols(1, data.n_cols - 1);
    arma::vec y = data.col(0);
    Environment fastglm = Environment::namespace_env("fastglm");
    Function fastglm_ = fastglm["fastglm"];
    List out;
    if (start.isNull()) {
      out = fastglm_(x, y, family);
    } else {
      arma::colvec start_ = as<arma::colvec>(start);
      out = fastglm_(x, y, family, Named("start") = start_);
    }
    arma::vec par = as<arma::vec>(out["coefficients"]);
    arma::vec residuals = as<arma::vec>(out["residuals"]);
    double value = out["deviance"];
    return List::create(Named("par") = par,
                  Named("value") = value / 2,
                  Named("residuals") = residuals);
  } else if (family == "lasso" || family == "gaussian") {
    arma::colvec theta_nonnull = as<arma::colvec>(theta);
    // Calculate negative log likelihood in gaussian family
    arma::vec y = data.col(0);
    arma::mat x = data.cols(1, data.n_cols - 1);
    double penalty = lambda * arma::accu(arma::abs(theta_nonnull));
    return List::create(
      Named("value") =
          arma::accu(arma::pow(y - x * theta_nonnull, 2)) / 2 + penalty
    );
  } else if (family == "binomial") {
    // Calculate negative log likelihood in binomial family
    arma::colvec theta_nonnull = as<arma::colvec>(theta);
    arma::vec y = data.col(0);
    arma::mat x = data.cols(1, data.n_cols - 1);
    arma::colvec u = x * theta_nonnull;
    return List::create(
      Named("value") = arma::accu(-y % u + arma::log(1 + exp(u)))
    );
  } else {
    // Calculate negative log likelihood in poisson family
    arma::colvec theta_nonnull = as<arma::colvec>(theta);
    arma::vec y = data.col(0);
    arma::mat x = data.cols(1, data.n_cols - 1);
    arma::colvec u = x * theta_nonnull;

    arma::colvec y_factorial(y.n_elem);
    for (unsigned int i = 0; i < y.n_elem; i++) {
      double log_factorial = 0;
      for (int j = 1; j <= y(i); ++j) {
        log_factorial += std::log(j);
      }
      y_factorial(i) = log_factorial;
    }

    return List::create(
      Named("value") = arma::accu(-y % u + exp(u) + y_factorial)
    );
  }
}

arma::colvec cost_update_gradient(
    arma::mat data,
    arma::colvec theta,
    std::string family
) {
  arma::rowvec new_data = data.row(data.n_rows - 1);
  arma::rowvec x = new_data.tail(new_data.n_elem - 1);
  double y = new_data(0);
  arma::colvec gradient;
  if (family.compare("binomial") == 0) {
    gradient = - (y - 1 / (1 + exp(-arma::as_scalar(x * theta)))) * x.t();
  } else if (family.compare("poisson") == 0) {
    gradient = - (y - exp(arma::as_scalar(x * theta))) * x.t();
  } else {
    // `family` is either "lasso" or "gaussian".
    gradient = - (y - arma::as_scalar(x * theta)) * x.t();
  }
  return gradient;
}

arma::mat cost_update_hessian(
    arma::mat data,
    arma::colvec theta,
    std::string family,
    double min_prob
) {
  arma::rowvec new_data = data.row(data.n_rows - 1);
  arma::rowvec x = new_data.tail(new_data.n_elem - 1);
  arma::mat hessian;
  if (family.compare("binomial") == 0) {
    double prob = 1 / (1 + exp(-arma::as_scalar(x * theta)));
    hessian = (x.t() * x) * arma::as_scalar((1 - prob) * prob);
  } else if (family.compare("poisson") == 0) {
    double prob = exp(arma::as_scalar(x * theta));
    hessian = (x.t() * x) * std::min(arma::as_scalar(prob), min_prob);
  } else {
    // `family` is either "lasso" or "gaussian".
    hessian = x.t() * x;
  }
  return hessian;
}

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
    std::function<arma::colvec(
        arma::mat data,
        arma::colvec theta,
        std::string family
    )> cost_gradient_wrapper,
    std::function<arma::mat(
        arma::mat data,
        arma::colvec theta,
        std::string family,
        double min_prob
    )> cost_hessian_wrapper
) {
  // Get the hessian
  arma::mat hessian_i = hessian.slice(i - 1);
  arma::colvec gradient;

  if (CUSTOM_FAMILIES.find(family) != CUSTOM_FAMILIES.end()) {
    arma::mat cost_hessian_result = cost_hessian_wrapper(
      data, theta_hat.col(i - 1),
      family,  // UNUSED
      min_prob  // UNUSED
    );
    arma::colvec cost_gradient_result = cost_gradient_wrapper(
      data, theta_hat.col(i - 1),
      family  // UNUSED
    );
    hessian_i += cost_hessian_result;
    gradient = cost_gradient_result;
  } else {
    hessian_i += cost_update_hessian(
      data, theta_hat.col(i - 1), family, min_prob
    );
    gradient = cost_update_gradient(data, theta_hat.col(i - 1), family);
  }

  // Add epsilon to the diagonal for PSD hessian
  arma::mat hessian_psd = hessian_i + epsilon * arma::eye<arma::mat>(
    theta_hat.n_rows, theta_hat.n_rows
  );

  // Calculate momentum step
  arma::vec momentum_step = arma::solve(hessian_psd, gradient);
  momentum = momentum_coef * momentum - momentum_step;

  // Update theta_hat with momentum
  theta_hat.col(i - 1) += momentum;

  // Winsorize if family is Poisson
  if (family == "poisson") {
    Environment desc_tools = Environment::namespace_env("DescTools");
    Function winsorize = desc_tools["Winsorize"];
    NumericVector winsorize_result = winsorize(
      Rcpp::_["x"] = theta_hat.col(i - 1),
      Rcpp::_["minval"] = winsorise_minval,
      Rcpp::_["maxval"] = winsorise_maxval
    );
    theta_hat.col(i - 1) = arma::vec(
      winsorize_result.begin(), winsorize_result.size(), false
    );
  } else if (family == "lasso" || family == "gaussian") {
    // Update theta_hat with L1 penalty
    double hessian_norm = arma::norm(hessian_i, "fro");
    arma::vec normd = arma::abs(theta_hat.col(i - 1)) - lambda / hessian_norm;
    theta_hat.col(i - 1) =
        arma::sign(theta_hat.col(i - 1)) % arma::max(
          normd, arma::zeros<arma::colvec>(normd.n_elem)
        );
  }

  for (int kk = 1; kk <= as<int>(k(data.n_rows - tau)); kk++) {
    for (unsigned j = tau + 1; j <= data.n_rows; j++) {
      if (CUSTOM_FAMILIES.find(family) != CUSTOM_FAMILIES.end()) {
        arma::mat cost_hessian_result = cost_hessian_wrapper(
          data.rows(tau, j - 1), theta_hat.col(i - 1),
          family,  // UNUSED
          min_prob  // UNUSED
        );
        hessian_i += cost_hessian_result;
        arma::colvec cost_gradient_result = cost_gradient_wrapper(
          data.rows(tau, j - 1), theta_hat.col(i - 1),
          family  // UNUSED
        );
        gradient = cost_gradient_result;
      } else {
        hessian_i += cost_update_hessian(
          data.rows(tau, j - 1), theta_hat.col(i - 1), family, min_prob
        );
        gradient = cost_update_gradient(
          data.rows(tau, j - 1), theta_hat.col(i - 1), family
        );
      }

      hessian_psd = hessian_i + epsilon * arma::eye<arma::mat>(
            theta_hat.n_rows, theta_hat.n_rows
          );
      momentum_step = arma::solve(hessian_psd, gradient);
      momentum = momentum_coef * momentum - momentum_step;
      theta_hat.col(i - 1) += momentum;

      // Winsorize if family is Poisson
      if (family == "poisson") {
        Environment desc_tools = Environment::namespace_env("DescTools");
        Function winsorize = desc_tools["Winsorize"];
        NumericVector winsorize_result = winsorize(
          Rcpp::_["x"] = theta_hat.col(i - 1),
          Rcpp::_["minval"] = winsorise_minval,
          Rcpp::_["maxval"] = winsorise_maxval
        );
        theta_hat.col(i - 1) =
            arma::vec(winsorize_result.begin(), winsorize_result.size(), false);
      } else if (family == "lasso" || family == "gaussian") {
        double hessian_norm = arma::norm(hessian_i, "fro");
        arma::vec normd = arma::abs(
          theta_hat.col(i - 1)
        ) - lambda / hessian_norm;
        theta_hat.col(i - 1) = arma::sign(theta_hat.col(i - 1)) % arma::max(
          normd, arma::zeros<arma::colvec>(normd.n_elem)
        );
      }
    }
  }

  theta_sum.col(i - 1) += theta_hat.col(i - 1);
  hessian.slice(i - 1) = std::move(hessian_i);
  return List::create(
    theta_hat.col(i - 1), theta_sum.col(i - 1), hessian.slice(i - 1), momentum
  );
}

List cost_optim(
    const std::string family,
    const int p,
    const arma::mat data_segment,
    Function cost,
    const double lambda,
    const bool cv
) {
  if (family == "vanilla") {
    return List::create(
      Named("par") = R_NilValue,
      Named("value") = cost(data_segment),
      Named("residuals") = R_NilValue
    );
  } else if (p == 1) {
    Environment stats = Environment::namespace_env("stats");
    Function optim = stats["optim"];
    List optim_result = optim(
      Named("par") = 0,
      Named("fn") = InternalFunction(
        +[](double theta, arma::mat data, Function cost) {
          return cost(
            Named("data") = data,
            Named("theta") = log(theta / (1 - theta))
          );
        }
      ),
      Named("method") = "Brent",
      Named("lower") = 0,
      Named("upper") = 1,
      Named("data") = data_segment,
      Named("cost") = cost
    );
    return List::create(
      Named("par") = log(
        as<double>(optim_result["par"]) / (1 - as<double>(optim_result["par"]))
      ),
      Named("value") = exp(as<double>(optim_result["value"])) /
        (1 + exp(as<double>(optim_result["value"]))),
      Named("residuals") = R_NilValue
    );
  } else if (p > 1) {
    Environment stats = Environment::namespace_env("stats");
    Function optim = stats["optim"];
    List optim_result = optim(
      Named("par") = arma::zeros<arma::vec>(p),
      Named("fn") = cost,
      Named("method") = "L-BFGS-B",
      Named("data") = data_segment
    );
    return List::create(
      Named("par") = optim_result["par"],
      Named("value") = optim_result["value"],
      Named("residuals") = R_NilValue
    );
  } else {
    // NOTEST
    Rcpp::stop("Internal error, this branch should not be reached.");
  }
}

List fastcpd_impl(
    arma::mat data,
    double beta,
    const int segment_count,
    const double trim,
    const double momentum_coef,
    Function k,
    const std::string family,
    const double epsilon,
    const double min_prob,
    const double winsorise_minval,
    const double winsorise_maxval,
    const int p,
    Nullable<Function> cost,
    Nullable<Function> cost_gradient,
    Nullable<Function> cost_hessian,
    const bool cp_only,
    const double vanilla_percentage
) {
  // Set up the initial values.
  const int n = data.n_rows;
  double lambda = 0;

  // After t = 1, the r_t_set R_t contains 0 and 1.
  arma::ucolvec r_t_set = {0, 1};
  // C(0) = NULL, C(1) = {0}.
  std::vector<arma::colvec> cp_sets = {arma::zeros<arma::vec>(0)};
  arma::linspace(1, n, n).for_each([&](int i) {
    cp_sets.push_back(arma::zeros<arma::vec>(1));
  });
  // Objective function: F(0) = -beta.
  arma::colvec f_t = arma::zeros<arma::vec>(n + 1);
  f_t(0) = -beta;

  fastcpd::FastcpdParameters fastcpd_parameters_class(
    data, beta, p, family, segment_count,
    winsorise_minval, winsorise_maxval, epsilon
  );

  fastcpd_parameters_class.wrap_cost(cost);
  fastcpd_parameters_class.wrap_cost_gradient(cost_gradient);
  fastcpd_parameters_class.wrap_cost_hessian(cost_hessian);

  if (vanilla_percentage != 1) {
    fastcpd_parameters_class.create_segment_indices();
    fastcpd_parameters_class.create_segment_statistics();
    fastcpd_parameters_class.update_beta();
    fastcpd_parameters_class.create_gradients();
  }

  for (int t = 2; t <= n; t++) {
    unsigned int r_t_count = r_t_set.n_elem;

    // Number of cost values is the same as the number of elements in R_t.
    arma::colvec cval = arma::zeros<arma::vec>(r_t_count);

    // For tau in R_t \ {t-1}.
    for (unsigned int i = 1; i < r_t_count; i++) {
      int tau = r_t_set(i - 1);
      if (family == "lasso") {
        // Mean of `err_sd` only works if error sd is unchanged.
        lambda = mean(
          fastcpd_parameters_class.get_err_sd()
        ) * sqrt(2 * log(p) / (t - tau));
      }
      arma::mat data_segment = data.rows(tau, t - 1);
      if (vanilla_percentage == 1 || t <= vanilla_percentage * n) {
        List cost_optim_result;
        if (CUSTOM_FAMILIES.find(family) != CUSTOM_FAMILIES.end()) {
          cost_optim_result = cost_optim(
            family, p, data_segment,
            fastcpd_parameters_class.cost.get(), 0, true
          );
        } else {
          cost_optim_result = negative_log_likelihood(
              data_segment, R_NilValue, family, 0, true
          );
        }
        cval(i - 1) = as<double>(cost_optim_result["value"]);
        if (vanilla_percentage < 1 && t <= vanilla_percentage * n) {
          fastcpd_parameters_class.update_theta_hat(
            i - 1, as<arma::colvec>(cost_optim_result["par"])
          );
          fastcpd_parameters_class.update_theta_sum(
            i - 1, as<arma::colvec>(cost_optim_result["par"])
          );
        }
      } else {
        List cost_update_result = cost_update(
          data.rows(0, t - 1),
          fastcpd_parameters_class.get_theta_hat(),
          fastcpd_parameters_class.get_theta_sum(),
          fastcpd_parameters_class.get_hessian(),
          tau,
          i,
          k,
          family,
          fastcpd_parameters_class.get_momentum(),
          momentum_coef,
          epsilon,
          min_prob,
          winsorise_minval,
          winsorise_maxval,
          lambda,
          fastcpd_parameters_class.cost_gradient_wrapper,
          fastcpd_parameters_class.cost_hessian_wrapper
        );
        fastcpd_parameters_class.update_theta_hat(
          i - 1, as<arma::colvec>(cost_update_result[0])
        );
        fastcpd_parameters_class.create_theta_sum(
          i - 1, as<arma::colvec>(cost_update_result[1])
        );
        fastcpd_parameters_class.update_hessian(
          i - 1, as<arma::mat>(cost_update_result[2])
        );
        fastcpd_parameters_class.update_momentum(
          as<arma::colvec>(cost_update_result[3])
        );
        arma::colvec theta =
            fastcpd_parameters_class.get_theta_sum().col(i - 1) / (t - tau);
        if (family == "poisson" && t - tau >= p) {
          Environment desc_tools = Environment::namespace_env("DescTools");
          Function winsorize = desc_tools["Winsorize"];
          NumericVector winsorize_result = winsorize(
            Rcpp::_["x"] = theta,
            Rcpp::_["minval"] = winsorise_minval,
            Rcpp::_["maxval"] = winsorise_maxval
          );
          theta = as<arma::colvec>(winsorize_result);
        }
        if (
          (
            FASTCPD_FAMILIES.find(family) != FASTCPD_FAMILIES.end() &&
            family != "lasso" && t - tau >= p
          ) ||
          (family == "lasso" && t - tau >= 3)
        ) {
          List cost_result = negative_log_likelihood(
            data_segment, Rcpp::wrap(theta), family, lambda
          );
          cval(i - 1) = as<double>(cost_result["value"]);
        } else if (CUSTOM_FAMILIES.find(family) != CUSTOM_FAMILIES.end()) {
          // if (warm_start && t - tau >= 50) {
          //   cost_result <- cost(data_segment, start = start[, tau + 1])
          //   start[, tau + 1] <- cost_result$par
          //   cval[i] <- cost_result$value
          // } else {
          Function cost_non_null = fastcpd_parameters_class.cost.get();
          SEXP cost_result = cost_non_null(data_segment, theta);
          cval(i - 1) = as<double>(cost_result);
          // }
        }
      }
    }

    if (vanilla_percentage != 1) {
      fastcpd_parameters_class.update_fastcpd_parameters(t);
    }

    // Step 3
    cval(r_t_count - 1) = 0;

    // `beta` adjustment seems to work but there might be better choices.
    arma::colvec obj = cval + f_t.rows(r_t_set) + beta;
    double min_obj = arma::min(obj);
    double tau_star = r_t_set(arma::index_min(obj));

    // Step 4
    cp_sets[t] = arma::join_cols(cp_sets[tau_star], arma::colvec{tau_star});

    // Step 5
    arma::ucolvec pruned_left =
        arma::find(cval + f_t.rows(r_t_set) <= min_obj);
    arma::ucolvec pruned_r_t_set =
        arma::zeros<arma::ucolvec>(pruned_left.n_elem + 1);
    pruned_r_t_set.rows(0, pruned_left.n_elem - 1) = r_t_set(pruned_left);
    pruned_r_t_set(pruned_left.n_elem) = t;
    r_t_set = std::move(pruned_r_t_set);

    if (vanilla_percentage != 1) {
      fastcpd_parameters_class.update_theta_hat(pruned_left);
      fastcpd_parameters_class.update_theta_sum(pruned_left);
      fastcpd_parameters_class.update_hessian(pruned_left);
    }

    // Objective function F(t).
    f_t(t) = min_obj;
  }

  // Remove change points close to the boundaries.
  arma::colvec cp_set = cp_sets[n];
  cp_set = cp_set(arma::find(cp_set >= trim * n));
  cp_set = cp_set(arma::find(cp_set <= (1 - trim) * n));
  arma::colvec cp_set_ = arma::zeros<arma::vec>(cp_set.n_elem + 1);
  if (cp_set.n_elem) {
    cp_set_.rows(1, cp_set_.n_elem - 1) = std::move(cp_set);
  }
  cp_set = arma::sort(arma::unique(std::move(cp_set_)));

  // Remove change points close to each other.
  arma::ucolvec cp_set_too_close = arma::find(arma::diff(cp_set) <= trim * n);
  if (cp_set_too_close.n_elem > 0) {
    int rest_element_count = cp_set.n_elem - cp_set_too_close.n_elem;
    arma::colvec cp_set_rest_left = arma::zeros<arma::vec>(rest_element_count),
                cp_set_rest_right = arma::zeros<arma::vec>(rest_element_count);
    for (unsigned int i = 0, i_left = 0, i_right = 0; i < cp_set.n_elem; i++) {
      if (
        arma::ucolvec left_find = arma::find(cp_set_too_close == i);
        left_find.n_elem == 0
      ) {
        cp_set_rest_left(i_left) = cp_set(i);
        i_left++;
      }
      if (
        arma::ucolvec right_find = arma::find(cp_set_too_close == i - 1);
        right_find.n_elem == 0
      ) {
        cp_set_rest_right(i_right) = cp_set(i);
        i_right++;
      }
    }
    cp_set = arma::floor((cp_set_rest_left + cp_set_rest_right) / 2);
  }
  cp_set = cp_set(arma::find(cp_set > 0));

  if (cp_only) {
    return List::create(
      Named("cp_set") = cp_set,
      Named("cost_values") = R_NilValue,
      Named("residual") = R_NilValue,
      Named("thetas") = R_NilValue
    );
  }

  arma::colvec cp_loc_ = arma::zeros<arma::colvec>(cp_set.n_elem + 2);
  if (cp_set.n_elem) {
    cp_loc_.rows(1, cp_loc_.n_elem - 2) = cp_set;
  }
  cp_loc_(cp_loc_.n_elem - 1) = n;
  arma::colvec cp_loc = arma::unique(std::move(cp_loc_));
  arma::colvec cost_values = arma::zeros<arma::vec>(cp_loc.n_elem - 1);
  arma::mat thetas = arma::zeros<arma::mat>(p, cp_loc.n_elem - 1);
  arma::colvec residual = arma::zeros<arma::colvec>(n);
  unsigned int residual_next_start = 0;
  for (unsigned int i = 0; i < cp_loc.n_elem - 1; i++) {
    arma::colvec segment_data_index_ =
        arma::linspace(cp_loc(i), cp_loc(i + 1) - 1, cp_loc(i + 1) - cp_loc(i));
    arma::ucolvec segment_data_index =
        arma::conv_to<arma::ucolvec>::from(std::move(segment_data_index_));
    arma::mat data_segment = data.rows(segment_data_index);
    List cost_optim_result;
    if (CUSTOM_FAMILIES.find(family) != CUSTOM_FAMILIES.end()) {
      cost_optim_result = cost_optim(
        family, p, data_segment,
        fastcpd_parameters_class.cost.get(), lambda, false
      );
    } else {
      cost_optim_result = negative_log_likelihood(
        data_segment, R_NilValue, family, lambda, false
      );
    }
    arma::colvec cost_optim_par = as<arma::colvec>(cost_optim_result["par"]);
    double cost_optim_value = as<double>(cost_optim_result["value"]);
    arma::colvec cost_optim_residual =
        as<arma::colvec>(cost_optim_result["residuals"]);
    thetas.col(i) = cost_optim_par;
    cost_values(i) = cost_optim_value;
    residual.rows(
      residual_next_start, residual_next_start + cost_optim_residual.n_elem - 1
    ) = cost_optim_residual;
    residual_next_start += cost_optim_residual.n_elem;
  }
  return List::create(
    Named("cp_set") = cp_set,
    Named("cost_values") = cost_values,
    Named("residual") = residual,
    Named("thetas") = thetas
  );
}
