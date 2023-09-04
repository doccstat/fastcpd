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
  hessian.slice(i - 1) = std::move(hessian_i);
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
    Function cost_gradient,
    Function cost_hessian,
    arma::ucolvec r_t_set,
    const int p,
    const double momentum_coef,
    const double min_prob,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double epsilon
) {
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
  arma::mat theta_hat(p, 1), theta_sum(p, 1);
  arma::cube hessian;
  if (family == "binomial") {
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    double prob(1 / (1 + exp(
      -arma::as_scalar(theta_hat.t() * data.row(0).tail(data.n_cols - 1).t())
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
    theta_hat.col(0) = arma::vec(
      winsorize_result.begin(), winsorize_result.size(), false
    );
    theta_sum.col(0) = arma::vec(
      winsorize_result.begin(), winsorize_result.size(), false
    );
    hessian = arma::cube(p, p, 1);
    hessian.slice(0) = (
      data.row(0).tail(data.n_cols - 1).t() * data.row(0).tail(data.n_cols - 1)
    ) * arma::as_scalar(
      exp(theta_hat.t() * data.row(0).tail(data.n_cols - 1).t())
    );
  } else if (family == "lasso" || family == "gaussian") {
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    hessian = arma::cube(p, p, 1);
    hessian.slice(0) = epsilon * arma::eye<arma::mat>(p, p) +
      data.row(0).tail(data.n_cols - 1).t() * data.row(0).tail(data.n_cols - 1);
  } else if (family == "custom") {
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    hessian = arma::cube(p, p, 1);
    hessian.slice(0) = arma::zeros<arma::mat>(p, p);
  }

  return List::create(
      Named("theta_hat") = theta_hat,
      Named("theta_sum") = theta_sum,
      Named("hessian") = hessian
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
  Environment base = Environment::namespace_env("base");
  Function formals = base["formals"],
            length = base["length"];
  if (as<unsigned int>(length(formals(cost))) == 1) {
    return List::create(
      Named("par") = R_NilValue,
      Named("value") = cost(data_segment),
      Named("residuals") = R_NilValue
    );
  } else if (family == "custom" && p == 1) {
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
  } else if (family == "custom" && p > 1) {
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
    return cost(data_segment, R_NilValue, family, lambda, cv);
  }
}

List init_fastcpd_parameters(
    const arma::mat data,
    const int p,
    const std::string family,
    const int segment_count,
    Function cost,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double epsilon,
    const double vanilla_percentage,
    double& beta
) {
  // `error_sd` is used in Gaussian family only. `act_num` is used in Lasso
  // and Gaussian families only.
  List fastcpd_parameters = List::create(
    Named("segment_indices") = R_NilValue,
    Named("segment_theta_hat") = R_NilValue,
    Named("err_sd") = arma::zeros<arma::vec>(segment_count),
    Named("act_num") = arma::zeros<arma::vec>(segment_count),
    Named("theta_hat") = arma::zeros<arma::mat>(p, 1),
    Named("theta_sum") = arma::zeros<arma::mat>(p, 1),
    Named("hessian") = R_NilValue,
    // Momentum will be used in the update step if `momentum_coef` is not 0.
    Named("momentum") = arma::zeros<arma::vec>(p)
  );
  if (vanilla_percentage != 1) {
    const int n = data.n_rows;

    // Segment the whole data set evenly based on the number of segments
    // specified in the `segment_count` parameter.
    unsigned int segment_length = floor(n / segment_count);
    int segment_remainder = n % segment_count;
    arma::colvec segment_indices = arma::zeros<arma::vec>(n);
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
    fastcpd_parameters["segment_indices"] = segment_indices;

    // Create a matrix to store the estimated coefficients in each segment,
    // where each row represents estimated coefficients for a segment.
    arma::mat segment_theta_hat = arma::zeros<arma::mat>(segment_count, p);
    arma::colvec err_sd = arma::zeros<arma::vec>(segment_count),
                act_num = arma::zeros<arma::vec>(segment_count);

    // Initialize theta_hat_t_t to be the estimate in the segment.
    for (
      int segment_index = 0; segment_index < segment_count; ++segment_index
    ) {
      arma::ucolvec segment_indices_ =
          arma::find(segment_indices == segment_index + 1);
      arma::mat data_segment = data.rows(segment_indices_);
      arma::rowvec segment_theta = as<arma::rowvec>(
        cost_optim(family, p, data_segment, cost, 0, false)["par"]
      );

      // Initialize the estimated coefficients for each segment to be the
      // estimated coefficients in the segment.
      segment_theta_hat.row(segment_index) = segment_theta;
      if (family == "lasso" || family == "gaussian") {
        arma::colvec segment_residual = data_segment.col(0) -
          data_segment.cols(1, data_segment.n_cols - 1) * segment_theta.t();
          double err_var =
              arma::as_scalar(arma::mean(arma::pow(segment_residual, 2)));
          err_sd(segment_index) = sqrt(err_var);
          act_num(segment_index) = arma::accu(arma::abs(segment_theta) > 0);
      }
    }
    fastcpd_parameters["segment_theta_hat"] = segment_theta_hat;
    fastcpd_parameters["err_sd"] = err_sd;
    fastcpd_parameters["act_num"] = act_num;

    // Adjust `beta` for Lasso and Gaussian families. This seems to be working
    // but there might be better choices.
    if (family == "lasso" || family == "gaussian") {
      beta = beta * (1 + mean(act_num));
    }

    List theta_hat_sum_hessian = init_theta_hat_sum_hessian(
      family, segment_theta_hat, data, p, winsorise_minval, winsorise_maxval,
      epsilon
    );
    fastcpd_parameters["theta_hat"] = theta_hat_sum_hessian["theta_hat"];
    fastcpd_parameters["theta_sum"] = theta_hat_sum_hessian["theta_sum"];
    fastcpd_parameters["hessian"] = theta_hat_sum_hessian["hessian"];
  }

  return fastcpd_parameters;
}

List append_fastcpd_parameters(
    List fastcpd_parameters,
    const double vanilla_percentage,
    const arma::mat data,
    const int t,
    const std::string family,
    const double winsorise_minval,
    const double winsorise_maxval,
    const int p,
    const double epsilon
) {
  if (vanilla_percentage != 1) {
    // for tau = t-1
    arma::rowvec new_data = data.row(t - 1).tail(data.n_cols - 1);
    arma::vec segment_indices = fastcpd_parameters["segment_indices"];
    const int segment_index = segment_indices(t - 1);
    arma::mat segment_theta_hat = fastcpd_parameters["segment_theta_hat"];
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
    } else if (family == "custom") {
      hessian_new = arma::zeros<arma::mat>(p, p);
    }

    arma::mat theta_hat = fastcpd_parameters["theta_hat"],
              theta_sum = fastcpd_parameters["theta_sum"];
    arma::cube hessian = fastcpd_parameters["hessian"];
    fastcpd_parameters["theta_hat"] = arma::join_rows(theta_hat, coef_add.t());
    fastcpd_parameters["theta_sum"] =
        arma::join_rows(theta_sum, cum_coef_add.t());
    fastcpd_parameters["hessian"] = arma::join_slices(hessian, hessian_new);
  }
  return fastcpd_parameters;
}
