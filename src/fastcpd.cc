#include "fastcpd.h"

using ::Rcpp::as;
using ::Rcpp::Environment;
using ::Rcpp::Function;
using ::Rcpp::List;
using ::Rcpp::Named;
using ::Rcpp::Nullable;
using ::Rcpp::NumericMatrix;
using ::Rcpp::NumericVector;
using ::Rcpp::S4;

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
    Function cost_gradient,
    Function cost_hessian
) {
  // Get the hessian
  arma::mat hessian_i = hessian.slice(i - 1);
  arma::colvec gradient;

  if (family == "custom") {
    NumericMatrix cost_hessian_result = cost_hessian(
      data, theta_hat.col(i - 1)
    );
    NumericVector cost_gradient_result = cost_gradient(
      data, theta_hat.col(i - 1)
    );
    hessian_i += arma::mat(
      cost_hessian_result.begin(),
      cost_hessian_result.nrow(),
      cost_hessian_result.ncol()
    );
    gradient = arma::colvec(
      cost_gradient_result.begin(), cost_gradient_result.size(), false
    );
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
      if (family == "custom") {
        NumericMatrix cost_hessian_result = cost_hessian(
          data.rows(tau, j - 1), theta_hat.col(i - 1)
        );
        hessian_i += arma::mat(
          cost_hessian_result.begin(),
          cost_hessian_result.nrow(),
          cost_hessian_result.ncol()
        );
        NumericVector cost_gradient_result = cost_gradient(
          data.rows(tau, j - 1), theta_hat.col(i - 1)
        );
        gradient = arma::colvec(
          cost_gradient_result.begin(), cost_gradient_result.size(), false
        );
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
  hessian.slice(i - 1) = hessian_i;
  return List::create(
    theta_hat.col(i - 1), theta_sum.col(i - 1), hessian.slice(i - 1), momentum
  );
}

List update_fastcpd_parameters(
    List fastcpd_parameters,
    arma::mat data,
    const int t,
    const int i,
    Function k,
    const int tau,
    const double lambda,
    const std::string family,
    const double vanilla_percentage,
    Function cost_gradient,
    Function cost_hessian,
    arma::vec r_t_set,
    const int p,
    const double momentum_coef,
    const double min_prob,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double epsilon
) {
  if (vanilla_percentage != 1) {
    List cost_update_result = cost_update(
      data.rows(0, t - 1),
      as<arma::mat>(fastcpd_parameters["theta_hat"]),
      as<arma::mat>(fastcpd_parameters["theta_sum"]),
      as<arma::cube>(fastcpd_parameters["hessian"]),
      tau,
      i,
      k,
      family,
      as<arma::colvec>(fastcpd_parameters["momentum"]),
      momentum_coef,
      epsilon,
      min_prob,
      winsorise_minval,
      winsorise_maxval,
      lambda,
      cost_gradient,
      cost_hessian
    );
    arma::mat theta_hat = fastcpd_parameters["theta_hat"],
              theta_sum = fastcpd_parameters["theta_sum"];
    arma::cube hessian = fastcpd_parameters["hessian"];
    theta_hat.col(i - 1) = as<arma::colvec>(cost_update_result[0]);
    theta_sum.col(i - 1) = as<arma::colvec>(cost_update_result[1]);
    hessian.slice(i - 1) = as<arma::mat>(cost_update_result[2]);
    fastcpd_parameters["theta_hat"] = theta_hat;
    fastcpd_parameters["theta_sum"] = theta_sum;
    fastcpd_parameters["hessian"] = hessian;
    fastcpd_parameters["momentum"] = cost_update_result[3];
  }
  return fastcpd_parameters;
}

List init_theta_hat_sum_hessian(
    const std::string family,
    const arma::mat segment_theta_hat,
    const arma::mat data,
    const int p,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double epsilon
) {
  arma::colvec theta_hat, theta_sum;
  arma::cube hessian;
  if (family == "binomial") {
    theta_sum = segment_theta_hat.row(0).t();
    theta_hat = segment_theta_hat.row(0).t();
    double prob = 1 / (1 + exp(-arma::as_scalar(
      theta_hat.t() * data.row(0).tail(data.n_cols - 1).t()
    )));
    hessian = arma::cube(p, p, 1);
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
    theta_hat = arma::vec(
      winsorize_result.begin(), winsorize_result.size(), false
    );
    theta_sum = arma::vec(
      winsorize_result.begin(), winsorize_result.size(), false
    );
    hessian = arma::cube(p, p, 1);
    hessian.slice(0) = (
      data.row(0).tail(data.n_cols - 1).t() * data.row(0).tail(data.n_cols - 1)
    ) * arma::as_scalar(
      exp(theta_hat.t() * data.row(0).tail(data.n_cols - 1).t())
    );
  } else if (family == "lasso" || family == "gaussian") {
    theta_hat = segment_theta_hat.row(0).t();
    theta_sum = segment_theta_hat.row(0).t();
    hessian = arma::cube(p, p, 1);
    hessian.slice(0) = epsilon * arma::eye<arma::mat>(p, p) +
      data.row(0).tail(data.n_cols - 1).t() * data.row(0).tail(data.n_cols - 1);
  } else if (family == "custom") {
    theta_hat = segment_theta_hat.row(0).t();
    theta_sum = segment_theta_hat.row(0).t();
    hessian = arma::cube(p, p, 1);
    hessian.slice(0) = arma::zeros<arma::mat>(p, p);
  }

  return List::create(
      Named("theta_hat") = theta_hat,
      Named("theta_sum") = theta_sum,
      Named("hessian") = hessian
  );
}
