#include "fastcpd_classes.h"
#include "fastcpd_constants.h"

namespace fastcpd::classes {

Fastcpd::Fastcpd(
    mat data,
    const double beta,
    const int p,
    const colvec order,
    const string family,
    const double vanilla_percentage,
    const int segment_count,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double epsilon,
    const double min_prob,
    const double momentum_coef,
    const colvec lower,
    const colvec upper,
    const mat mean_data_cov,
    const unsigned int p_response
) : data(data),
    beta(beta),
    p(p),
    order(order),
    family(family),
    vanilla_percentage(vanilla_percentage),
    segment_count(segment_count),
    winsorise_minval(winsorise_minval),
    winsorise_maxval(winsorise_maxval),
    epsilon(epsilon),
    min_prob(min_prob),
    momentum_coef(momentum_coef),
    lower(lower),
    upper(upper),
    mean_data_cov(mean_data_cov),
    p_response(p_response) {
  n = data.n_rows;
  segment_indices = vec(n);
  segment_theta_hat = mat(segment_count, p);
  err_sd = vec(segment_count);
  act_num = vec(segment_count);
  theta_hat = mat(p, 1);
  theta_sum = mat(p, 1);
  hessian = cube(p, p, 1);
  momentum = vec(p);

  if (family == "variance") {
    variance_data_mean = mean(data, 0);
  }

  // TODO(doccstat): Store environment functions from R.
}

colvec Fastcpd::get_err_sd() {
  return err_sd;
}

void Fastcpd::update_err_sd(
  const unsigned int segment_index, const double err_var
) {
  err_sd(segment_index) = sqrt(err_var);
}

void Fastcpd::create_segment_indices() {
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

void Fastcpd::update_theta_hat(
  const unsigned int col, colvec new_theta_hat
) {
  theta_hat.col(col) = new_theta_hat;
}

void Fastcpd::update_theta_hat(colvec new_theta_hat) {
  theta_hat = arma::join_rows(theta_hat, new_theta_hat);
}

void Fastcpd::create_theta_sum(
  const unsigned int col, colvec new_theta_sum
) {
  theta_sum.col(col) = new_theta_sum;
}

mat Fastcpd::get_theta_sum() {
  return theta_sum;
}

void Fastcpd::update_theta_sum(
  const unsigned int col, colvec new_theta_sum
) {
  theta_sum.col(col) += new_theta_sum;
}

void Fastcpd::update_theta_sum(colvec new_theta_sum) {
  theta_sum = arma::join_rows(theta_sum, new_theta_sum);
}

void Fastcpd::update_hessian(
  const unsigned int slice, mat new_hessian
) {
  hessian.slice(slice) = new_hessian;
}

void Fastcpd::update_hessian(mat new_hessian) {
  hessian = arma::join_slices(hessian, new_hessian);
}

void Fastcpd::update_theta_hat(ucolvec pruned_left) {
  theta_hat = theta_hat.cols(pruned_left);
}

void Fastcpd::update_theta_sum(ucolvec pruned_left) {
  theta_sum = theta_sum.cols(pruned_left);
}

void Fastcpd::update_hessian(ucolvec pruned_left) {
  hessian = hessian.slices(pruned_left);
}

// TODO(doccstat): Use `segment_theta` as warm start.

void Fastcpd::create_segment_statistics() {
  for (
    int segment_index = 0; segment_index < segment_count; ++segment_index
  ) {
    ucolvec segment_indices_ =
        arma::find(segment_indices == segment_index + 1);
    mat data_segment = data.rows(segment_indices_);
    rowvec segment_theta;
    if (!contain(FASTCPD_FAMILIES, family)) {
      segment_theta = as<rowvec>(cost_optim(
        data_segment, 0, false
      )["par"]);
    } else {
      segment_theta = as<rowvec>(
        cost_function_wrapper(
          data_segment, R_NilValue, 0, true, R_NilValue
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
  if (family == "lasso") {
    beta = beta * (1 + mean(act_num));
  }
}

double Fastcpd::get_beta() {
  return beta;
}

void Fastcpd::update_momentum(colvec new_momentum) {
  momentum = new_momentum;
}

void Fastcpd::create_gradients() {
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
    theta_hat.col(0) = clamp(
      segment_theta_hat.row(0).t(), winsorise_minval, winsorise_maxval
    );
    theta_sum.col(0) = clamp(
      segment_theta_hat.row(0).t(), winsorise_minval, winsorise_maxval
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
  } else if (!contain(FASTCPD_FAMILIES, family)) {
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    hessian.slice(0) = zeros<mat>(p, p);
  }
}

void Fastcpd::update_fastcpd_parameters(const unsigned int t) {
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
    coef_add = clamp(coef_add, winsorise_minval, winsorise_maxval);
    cum_coef_add = clamp(coef_add, winsorise_minval, winsorise_maxval);
    hessian_new =
        (new_data.t() * new_data) * as_scalar(
          exp(coef_add * new_data.t())
        );
  } else if (family == "lasso" || family == "gaussian") {
    hessian_new = new_data.t() * new_data + epsilon * eye<mat>(p, p);
  } else if (family == "arma") {
    hessian_new = cost_update_hessian(data.rows(0, t - 1), coef_add.t());
  } else if (!contain(FASTCPD_FAMILIES, family)) {
    hessian_new = zeros<mat>(p, p);
  }

  update_theta_hat(coef_add.t());
  update_theta_sum(cum_coef_add.t());
  update_hessian(hessian_new);
}

void Fastcpd::wrap_cost(Nullable<Function> cost) {
  this->cost = cost;
  if (contain(FASTCPD_FAMILIES, family)) {
    cost_function_wrapper = std::bind(  // # nocov start
      &Fastcpd::negative_log_likelihood,  // # nocov end
      this,
      std::placeholders::_1,
      std::placeholders::_2,
      std::placeholders::_3,
      std::placeholders::_4,
      std::placeholders::_5
    );
  } else {
    fastcpd::classes::CostFunction costFunction(cost.get());
    cost_function_wrapper = costFunction;
  }
}

void Fastcpd::wrap_cost_gradient(Nullable<Function> cost_gradient) {
  this->cost_gradient = cost_gradient;
  if (contain(FASTCPD_FAMILIES, family)) {
    cost_gradient_wrapper = std::bind(  // # nocov start
      &Fastcpd::cost_update_gradient,  // # nocov end
      this,
      std::placeholders::_1,
      std::placeholders::_2
    );
  } else if (cost_gradient.isNotNull()) {
    fastcpd::classes::CostGradient costGradient(cost_gradient.get());
    cost_gradient_wrapper = costGradient;
  } else if (cost_gradient.isNull()) {
    // `cost_gradient` can be `NULL` in the case of vanilla PELT.
  } else {
    // # nocov start
    stop("This branch should not be reached at classes.cc: 290.");
    // # nocov end
  }
}

void Fastcpd::wrap_cost_hessian(Nullable<Function> cost_hessian) {
  this->cost_hessian = cost_hessian;
  if (contain(FASTCPD_FAMILIES, family)) {
    cost_hessian_wrapper = std::bind(  // # nocov start
      &Fastcpd::cost_update_hessian,  // # nocov end
      this,
      std::placeholders::_1,
      std::placeholders::_2
    );
  } else if (cost_hessian.isNotNull()) {
    fastcpd::classes::CostHessian costHessian(cost_hessian.get());
    cost_hessian_wrapper = costHessian;
  } else if (cost_hessian.isNull()) {
    // `cost_hessian` can be `NULL` in the case of vanilla PELT.
  } else {
    // # nocov start
    stop("This branch should not be reached at classes.cc: 304.");
    // # nocov end
  }
}

void Fastcpd::cost_update(
  const unsigned int t,
  const unsigned int tau,
  const unsigned int i,
  Function k,
  const double lambda,
  const colvec line_search
) {
  List cost_update_result = cost_update_steps(
    data.rows(0, t - 1), tau, i, k, momentum, lambda, line_search
  );
  update_theta_hat(i - 1, as<colvec>(cost_update_result[0]));
  create_theta_sum(i - 1, as<colvec>(cost_update_result[1]));
  update_hessian(i - 1, as<mat>(cost_update_result[2]));
  update_momentum(as<colvec>(cost_update_result[3]));
}

List Fastcpd::cost_update_steps(
    const mat data,
    const int tau,
    const int i,
    Function k,
    colvec momentum,
    const double lambda,
    colvec line_search
) {
  cost_update_step(data, i, 0, data.n_rows - 1, lambda, line_search);

  for (int kk = 1; kk <= as<int>(k(data.n_rows - tau)); kk++) {
    for (unsigned j = tau + 1; j <= data.n_rows; j++) {
      cost_update_step(data, i, tau, j - 1, lambda, line_search);
    }
  }

  theta_sum.col(i - 1) += theta_hat.col(i - 1);
  return List::create(
    theta_hat.col(i - 1), theta_sum.col(i - 1), hessian.slice(i - 1), momentum
  );
}

void Fastcpd::cost_update_step(
  const mat data,
  const int i,
  const int data_start,
  const int data_end,
  const double lambda,
  const colvec line_search
) {
  mat hessian_i = hessian.slice(i - 1);
  colvec gradient;

  if (!contain(FASTCPD_FAMILIES, family)) {
    mat cost_hessian_result = cost_hessian_wrapper(
      data.rows(data_start, data_end), theta_hat.col(i - 1)
    );
    hessian_i += cost_hessian_result;
    colvec cost_gradient_result = cost_gradient_wrapper(
      data.rows(data_start, data_end), theta_hat.col(i - 1)
    );
    gradient = cost_gradient_result;
  } else {
    hessian_i += cost_update_hessian(
      data.rows(data_start, data_end), theta_hat.col(i - 1)
    );
    gradient = cost_update_gradient(
      data.rows(data_start, data_end), theta_hat.col(i - 1)
    );
  }

  // Add epsilon to the diagonal for PSD hessian
  mat hessian_psd =
    hessian_i + epsilon * eye<mat>(theta_hat.n_rows, theta_hat.n_rows);

  // Calculate momentum step
  momentum = momentum_coef * momentum - solve(hessian_psd, gradient);

  double best_learning_rate = 1;
  colvec line_search_costs = zeros<colvec>(line_search.n_elem);

  // Line search
  if (line_search.n_elem > 1 || line_search[0] != 1) {
    for (
      unsigned int line_search_index = 0;
      line_search_index < line_search.n_elem;
      line_search_index++
    ) {
      line_search_costs[line_search_index] = cost_function_wrapper(
        data, Rcpp::wrap(max(min(
          theta_hat.col(i - 1) + line_search[line_search_index] * momentum,
          upper
        ), lower)), lambda, false, R_NilValue
      )["value"];
    }
  }
  best_learning_rate = line_search[line_search_costs.index_min()];

  // Update theta_hat with momentum
  theta_hat.col(i - 1) += best_learning_rate * momentum;

  theta_hat.col(i - 1) = min(theta_hat.col(i - 1), upper);
  theta_hat.col(i - 1) = max(theta_hat.col(i - 1), lower);

  // Winsorize if family is Poisson
  if (family == "poisson") {
    theta_hat.col(i - 1) = clamp(
      theta_hat.col(i - 1), winsorise_minval, winsorise_maxval
    );
  } else if (family == "lasso" || family == "gaussian") {
    // Update theta_hat with L1 penalty
    double hessian_norm = norm(hessian_i, "fro");
    vec normd = arma::abs(theta_hat.col(i - 1)) - lambda / hessian_norm;
    theta_hat.col(i - 1) =
      sign(theta_hat.col(i - 1)) % max(normd, zeros<colvec>(normd.n_elem));
  }

  hessian.slice(i - 1) = std::move(hessian_i);
}

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

List Fastcpd::cost_optim(
    const mat data_segment,
    const double lambda,
    const bool cv
) {
  Function cost_ = cost.get();
  List cost_optim_result;
  if (vanilla_percentage == 1) {
    cost_optim_result = List::create(
      Named("par") = R_NilValue,
      Named("value") = cost_(data_segment),
      Named("residuals") = R_NilValue
    );
  } else if (p == 1) {
    Environment stats = Environment::namespace_env("stats");
    Function optim = stats["optim"];
    List optim_result = optim(
      Named("par") = 0,
      Named("fn") = InternalFunction(
        +[](double theta, mat data, Function cost_) {
          return cost_(
            Named("data") = data,
            Named("theta") = log(theta / (1 - theta))
          );
        }
      ),
      Named("method") = "Brent",
      Named("lower") = 0,
      Named("upper") = 1,
      Named("data") = data_segment,
      Named("cost") = cost_
    );
    cost_optim_result = List::create(
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
      Named("par") = zeros<vec>(p),
      Named("fn") = cost_,
      Named("method") = "L-BFGS-B",
      Named("data") = data_segment,
      Named("lower") = lower,
      Named("upper") = upper
    );
    cost_optim_result = List::create(
      Named("par") = optim_result["par"],
      Named("value") = optim_result["value"],
      Named("residuals") = R_NilValue
    );
  } else {
    // # nocov start
    stop("This branch should not be reached at classes.cc: 945.");
    // # nocov end
  }
  return cost_optim_result;
}

Nullable<colvec> Fastcpd::get_segment_theta_hat(const unsigned int t) {
  const unsigned int segment_index = segment_indices(t - 1);
  return Rcpp::wrap(segment_theta_hat[segment_index - 1]);
}

}  // namespace fastcpd::classes
