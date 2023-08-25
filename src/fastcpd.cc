#include "fastcpd.h"

Rcpp::List negative_log_likelihood(
    arma::mat data,
    Rcpp::Nullable<arma::colvec> theta,
    std::string family,
    double lambda,
    bool cv,
    Rcpp::Nullable<arma::colvec> start
) {
  if (theta.isNull() && family == "lasso" && cv) {
    Rcpp::Environment glmnet = Rcpp::Environment::namespace_env("glmnet");
    Rcpp::Function cv_glmnet = glmnet["cv.glmnet"], predict_glmnet = glmnet["predict.glmnet"];
    Rcpp::List out = cv_glmnet(
      data.cols(1, data.n_cols - 1),
      data.col(0),
      Rcpp::Named("family") = "gaussian"
    );
    Rcpp::S4 out_coef = predict_glmnet(
      out["glmnet.fit"],
      Rcpp::Named("s") = out["lambda.1se"],
      Rcpp::Named("type") = "coefficients",
      Rcpp::Named("exact") = false
    );
    arma::vec glmnet_i = Rcpp::as<arma::vec>(out_coef.slot("i"));
    arma::vec glmnet_x = Rcpp::as<arma::vec>(out_coef.slot("x"));
    arma::vec par = arma::zeros(data.n_cols - 1);
    for (unsigned int i = 1; i < glmnet_i.n_elem; i++) {
      par(glmnet_i(i) - 1) = glmnet_x(i);
    }
    return Rcpp::List::create(Rcpp::Named("par") = par,
                  Rcpp::Named("value") = R_NilValue);
  } else if (theta.isNull() && family == "lasso" && !cv) {
    Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats"), glmnet = Rcpp::Environment::namespace_env("glmnet");
    Rcpp::Function deviance = stats["deviance"], glmnet_ = glmnet["glmnet"], predict_glmnet = glmnet["predict.glmnet"];
    Rcpp::List out = glmnet_(
      data.cols(1, data.n_cols - 1),
      data.col(0),
      Rcpp::Named("family") = "gaussian",
      Rcpp::Named("lambda") = lambda
    );
    Rcpp::S4 out_par = out["beta"];
    arma::vec par_i = Rcpp::as<arma::vec>(out_par.slot("i"));
    arma::vec par_x = Rcpp::as<arma::vec>(out_par.slot("x"));
    arma::vec par = arma::zeros(data.n_cols - 1);
    for (unsigned int i = 0; i < par_i.n_elem; i++) {
      par(par_i(i)) = par_x(i);
    }
    double value = Rcpp::as<double>(deviance(out));
    arma::vec fitted_values = Rcpp::as<arma::vec>(predict_glmnet(out, data.cols(1, data.n_cols - 1), Rcpp::Named("s") = lambda));
    arma::vec residuals = data.col(0) - fitted_values;
    return Rcpp::List::create(Rcpp::Named("par") = par,
                  Rcpp::Named("value") = value / 2,
                  Rcpp::Named("residuals") = residuals);
  } else if (theta.isNull()) {
    // Estimate theta in binomial/poisson/gaussian family
    arma::mat x = data.cols(1, data.n_cols - 1);
    arma::vec y = data.col(0);
    Rcpp::Environment fastglm = Rcpp::Environment::namespace_env("fastglm");
    Rcpp::Function fastglm_ = fastglm["fastglm"];
    Rcpp::List out;
    if (start.isNull()) {
      out = fastglm_(x, y, family);
    } else {
      arma::colvec start_ = Rcpp::as<arma::colvec>(start);
      out = fastglm_(x, y, family, Rcpp::Named("start") = start_);
    }
    arma::vec par = Rcpp::as<arma::vec>(out["coefficients"]);
    arma::vec residuals = Rcpp::as<arma::vec>(out["residuals"]);
    double value = out["deviance"];
    return Rcpp::List::create(Rcpp::Named("par") = par,
                  Rcpp::Named("value") = value / 2,
                  Rcpp::Named("residuals") = residuals);
  } else if (family == "lasso" || family == "gaussian") {

    arma::colvec theta_nonnull = Rcpp::as<arma::colvec>(theta);
    // Calculate negative log likelihood in gaussian family
    arma::vec y = data.col(0);
    arma::mat x = data.cols(1, data.n_cols - 1);
    double penalty = lambda * arma::accu(arma::abs(theta_nonnull));
    return Rcpp::List::create(Rcpp::Named("value") = arma::accu(arma::pow(y - x * theta_nonnull, 2)) / 2 + penalty);

  } else if (family == "binomial") {

    // Calculate negative log likelihood in binomial family
    arma::colvec theta_nonnull = Rcpp::as<arma::colvec>(theta);
    arma::vec y = data.col(0);
    arma::mat x = data.cols(1, data.n_cols - 1);
    arma::colvec u = x * theta_nonnull;
    return Rcpp::List::create(Rcpp::Named("value") = arma::accu(-y % u + arma::log(1 + arma::exp(u))));
  } else {
    // Calculate negative log likelihood in poisson family
    arma::colvec theta_nonnull = Rcpp::as<arma::colvec>(theta);
    arma::vec y = data.col(0);
    arma::mat x = data.cols(1, data.n_cols - 1);
    arma::colvec u = x * theta_nonnull;
    Rcpp::NumericVector y_factorial = Rcpp::lfactorial(Rcpp::wrap(y));
    return Rcpp::List::create(Rcpp::Named("value") = arma::accu(-y % u + arma::exp(u) + arma::vec(y_factorial.begin(), y_factorial.size(), false)));
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

Rcpp::List cost_update(
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
    const double lambda,
    Rcpp::Function cost_gradient,
    Rcpp::Function cost_hessian
) {
  // Get the hessian
  arma::mat hessian_i = hessian.slice(i - 1);
  arma::colvec gradient;

  if (family == "custom") {
    Rcpp::NumericMatrix cost_hessian_result = cost_hessian(data, theta_hat.col(i - 1));
    Rcpp::NumericVector cost_gradient_result = cost_gradient(data, theta_hat.col(i - 1));
    hessian_i += arma::mat(cost_hessian_result.begin(), cost_hessian_result.nrow(), cost_hessian_result.ncol());
    gradient = arma::colvec(cost_gradient_result.begin(), cost_gradient_result.size(), false);
  } else {
    hessian_i += cost_update_hessian(data, theta_hat.col(i - 1), family, min_prob);
    gradient = cost_update_gradient(data, theta_hat.col(i - 1), family);
  }

  // Add epsilon to the diagonal for PSD hessian
  arma::mat hessian_psd = hessian_i + epsilon * arma::eye<arma::mat>(theta_hat.n_rows, theta_hat.n_rows);

  // Calculate momentum step
  arma::vec momentum_step = arma::solve(hessian_psd, gradient);
  momentum = momentum_coef * momentum - momentum_step;

  // Update theta_hat with momentum
  theta_hat.col(i - 1) += momentum;

  // Winsorize if family is Poisson
  if (family == "poisson") {
    Rcpp::Environment desc_tools = Rcpp::Environment::namespace_env("DescTools");
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
    theta_hat.col(i - 1) = arma::sign(theta_hat.col(i - 1)) % arma::max(normd, arma::zeros<arma::colvec>(normd.n_elem));
  }

  for (int kk = 1; kk <= Rcpp::as<int>(k(data.n_rows - tau)); kk++) {
    for (unsigned j = tau + 1; j <= data.n_rows; j++) {
      if (family == "custom") {
        Rcpp::NumericMatrix cost_hessian_result = cost_hessian(data.rows(tau, j - 1), theta_hat.col(i - 1));
        hessian_i += arma::mat(cost_hessian_result.begin(), cost_hessian_result.nrow(), cost_hessian_result.ncol());
        Rcpp::NumericVector cost_gradient_result = cost_gradient(data.rows(tau, j - 1), theta_hat.col(i - 1));
        gradient = arma::colvec(cost_gradient_result.begin(), cost_gradient_result.size(), false);
      } else {
        hessian_i += cost_update_hessian(data.rows(tau, j - 1), theta_hat.col(i - 1), family, min_prob);
        gradient = cost_update_gradient(data.rows(tau, j - 1), theta_hat.col(i - 1), family);
      }

      hessian_psd = hessian_i + epsilon * arma::eye<arma::mat>(theta_hat.n_rows, theta_hat.n_rows);
      momentum_step = arma::solve(hessian_psd, gradient);
      momentum = momentum_coef * momentum - momentum_step;
      theta_hat.col(i - 1) += momentum;

      // Winsorize if family is Poisson
      if (family == "poisson") {
        Rcpp::Environment desc_tools = Rcpp::Environment::namespace_env("DescTools");
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
        theta_hat.col(i - 1) = arma::sign(theta_hat.col(i - 1)) % arma::max(normd, arma::zeros<arma::colvec>(normd.n_elem));
      }
    }
  }

  theta_sum.col(i - 1) += theta_hat.col(i - 1);
  hessian.slice(i - 1) = hessian_i;
  return Rcpp::List::create(theta_hat.col(i - 1), theta_sum.col(i - 1), hessian.slice(i - 1), momentum);
}
