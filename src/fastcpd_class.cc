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
      segment_theta = as<rowvec>(cost_optim(data_segment, 0)["par"]);
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

List Fastcpd::cost_optim(
    const mat data_segment,
    const double lambda
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
