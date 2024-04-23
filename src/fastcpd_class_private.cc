#include "fastcpd_classes.h"
#include "fastcpd_constants.h"
#include "RProgress.h"

namespace fastcpd::classes {

void Fastcpd::create_clock_in_r(const std::string name) {
  if (!r_clock.empty()) {
    rClock.stop(name);
  }
}

void Fastcpd::create_cost_function_wrapper(Nullable<Function> cost) {
  DEBUG_RCOUT(family);
  if (contain(FASTCPD_FAMILIES, family)) {
    cost_function_wrapper = std::bind(  // # nocov start
      &Fastcpd::get_cost_result,  // # nocov end
      this,
      std::placeholders::_1,
      std::placeholders::_2,
      std::placeholders::_3,
      std::placeholders::_4,
      std::placeholders::_5,
      std::placeholders::_6
    );
  } else {
    fastcpd::classes::CostFunction costFunction(cost.get(), data);
    cost_function_wrapper = costFunction;
  }
}

void Fastcpd::create_cost_gradient_wrapper(Nullable<Function> cost_gradient) {
  if (contain(FASTCPD_FAMILIES, family)) {
    cost_gradient_wrapper = std::bind(  // # nocov start
      get_gradient,  // # nocov end
      this,
      std::placeholders::_1,
      std::placeholders::_2,
      std::placeholders::_3
    );
  } else if (cost_gradient.isNotNull()) {
    fastcpd::classes::CostGradient costGradient(cost_gradient.get(), data);
    cost_gradient_wrapper = costGradient;
  } else if (cost_gradient.isNull()) {
    // `cost_gradient` can be `NULL` in the case of vanilla PELT.
  } else {
    // # nocov start
    stop("This branch should not be reached at classes.cc: 290.");
    // # nocov end
  }
}

void Fastcpd::create_cost_hessian_wrapper(Nullable<Function> cost_hessian) {
  if (contain(FASTCPD_FAMILIES, family)) {
    cost_hessian_wrapper = std::bind(  // # nocov start
      get_hessian,  // # nocov end
      this,
      std::placeholders::_1,
      std::placeholders::_2,
      std::placeholders::_3
    );
  } else if (cost_hessian.isNotNull()) {
    fastcpd::classes::CostHessian costHessian(cost_hessian.get(), data);
    cost_hessian_wrapper = costHessian;
  } else if (cost_hessian.isNull()) {
    // `cost_hessian` can be `NULL` in the case of vanilla PELT.
  } else {
    // # nocov start
    stop("This branch should not be reached at classes.cc: 304.");
    // # nocov end
  }
}

void Fastcpd::create_gradients() {
  if (vanilla_percentage == 1) return;
  if (family == "binomial") {
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    const double prob = 1 / (1 + exp(-dot(theta_hat, data.row(0).tail(p))));
    hessian.slice(0) = (data.row(0).tail(p).t() * data.row(0).tail(p)) *
      prob * (1 - prob);
  } else if (family == "poisson") {
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    hessian.slice(0) = (
      data.row(0).tail(p).t() * data.row(0).tail(p)
    ) * exp(dot(theta_hat, data.row(0).tail(p)));
  } else if (family == "lasso" || family == "gaussian") {
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    hessian.slice(0) = data.row(0).tail(p).t() * data.row(0).tail(p) +
      epsilon * eye<mat>(p, p);
  } else if (!contain(FASTCPD_FAMILIES, family)) {
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    hessian.slice(0) = zeros<mat>(p, p);
  }
}

// TODO(doccstat): Use `segment_theta` as warm start.

void Fastcpd::create_segment_statistics() {
  if (!contain(FASTCPD_FAMILIES, family) && vanilla_percentage == 1) return;
  for (
    int segment_index = 0; segment_index < segment_count; ++segment_index
  ) {
    rowvec segment_theta;
    if (!contain(FASTCPD_FAMILIES, family)) {
      segment_theta = get_optimized_cost(
        segment_indices(segment_index), segment_indices(segment_index + 1) - 1
      ).par;
    } else {
      segment_theta = get_cost_result(
        segment_indices(segment_index),
        segment_indices(segment_index + 1) - 1,
        R_NilValue, 0,
        true,
        R_NilValue
      ).par;
    }

    // Initialize the estimated coefficients for each segment to be the
    // estimated coefficients in the segment.
    segment_theta_hat.row(segment_index) = segment_theta;
    if (family == "lasso" || family == "gaussian") {
      mat data_segment = data.rows(
        segment_indices(segment_index), segment_indices(segment_index + 1) - 1
      );
      colvec segment_residual = data_segment.col(0) -
        data_segment.cols(1, data_segment.n_cols - 1) * segment_theta.t();
      double err_var = as_scalar(mean(square(segment_residual)));
      update_err_sd(segment_index, err_var);
      act_num(segment_index) = accu(abs(segment_theta) > 0);
    }
  }
  if (family == "lasso") {
    beta = beta * (1 + mean(act_num));
  }
}

double Fastcpd::get_cost_adjustment_value(const unsigned nrows) {
  double adjusted = 0;
  if (cost_adjustment == "MBIC" || cost_adjustment == "MDL") {
    adjusted = p * std::log((double) nrows / data_n_rows) / 2.0;
  }
  if (cost_adjustment == "MDL") {
    adjusted *= std::log2(M_E);
  }
  return adjusted;
}

CostResult Fastcpd::get_cost_result(
    const unsigned int segment_start,
    const unsigned int segment_end,
    Nullable<colvec> theta,
    const double lambda,
    const bool cv,
    Nullable<colvec> start
) {
  CostResult cost_result;
  if (theta.isNull()) {
    cost_result = get_nll_wo_theta(
      segment_start, segment_end, lambda, cv, start
    );
  } else {
    cost_result = CostResult{
      {colvec()},
      {colvec()},
      get_nll_wo_cv(segment_start, segment_end, as<colvec>(theta), lambda)
    };
  }
  cost_result.value = update_cost_value(
    cost_result.value, segment_end - segment_start + 1
  );
  return cost_result;
}

List Fastcpd::get_cp_set(const colvec raw_cp_set, const double lambda) {
  colvec cp_set = update_cp_set(raw_cp_set);

  if (cp_only) {
    return List::create(
      Named("raw_cp_set") = raw_cp_set,
      Named("cp_set") = cp_set,
      Named("cost_values") = R_NilValue,
      Named("residual") = R_NilValue,
      Named("thetas") = R_NilValue
    );
  }

  colvec cp_loc_ = zeros<colvec>(cp_set.n_elem + 2);
  if (cp_set.n_elem) {
    cp_loc_.rows(1, cp_loc_.n_elem - 2) = cp_set;
  }
  cp_loc_(cp_loc_.n_elem - 1) = data_n_rows;
  colvec cp_loc = unique(std::move(cp_loc_));
  colvec cost_values = zeros<vec>(cp_loc.n_elem - 1);
  mat thetas = zeros<mat>(p, cp_loc.n_elem - 1);
  mat residual;
  if (
    family == "mean" || family == "variance" || family == "meanvariance"
  ) {
    residual = zeros<mat>(data_n_rows, d);
  } else if (family == "mgaussian") {
    residual = zeros<mat>(data_n_rows, p_response);
  } else {
    residual = zeros<mat>(data_n_rows, 1);
  }
  unsigned int residual_next_start = 0;

  for (unsigned int i = 0; i < cp_loc.n_elem - 1; i++) {
    CostResult cost_result;
    if (!contain(FASTCPD_FAMILIES, family)) {
      cost_result = get_optimized_cost(cp_loc(i), cp_loc(i + 1) - 1);
    } else {
      cost_result = get_cost_result(
        cp_loc(i), cp_loc(i + 1) - 1, R_NilValue, lambda, false, R_NilValue
      );
    }

    cost_values(i) = cost_result.value;

    // Parameters are not involved for PELT.
    if (vanilla_percentage < 1) {
      thetas.col(i) = colvec(cost_result.par);
    }

    // Residual is only calculated for built-in families.
    if (
      contain(FASTCPD_FAMILIES, family) &&
      !(family == "mean" || family == "variance" || family == "meanvariance")
    ) {
      mat cost_optim_residual = cost_result.residuals;
      residual.rows(
        residual_next_start,
        residual_next_start + cost_optim_residual.n_rows - 1
      ) = cost_optim_residual;
      residual_next_start += cost_optim_residual.n_rows;
    }
  }
  return List::create(
    Named("raw_cp_set") = raw_cp_set,
    Named("cp_set") = cp_set,
    Named("cost_values") = cost_values,
    Named("residual") = residual,
    Named("thetas") = thetas
  );
}

double Fastcpd::get_cval_for_r_t_set(
  const int tau,
  const unsigned int i,
  const int t,
  double lambda
) {
  DEBUG_RCOUT(i);
  if (family == "lasso") {
    // Mean of `err_sd` only works if error sd is unchanged.
    lambda = mean(err_sd) * sqrt(2 * std::log(p) / (t - tau));
  }
  if (t > vanilla_percentage * data_n_rows) {
    return get_cval_sen(tau, t - 1, i, lambda);
  } else {
    return get_cval_pelt(tau, t - 1, i, lambda);
  }
}

double Fastcpd::get_cval_pelt(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const unsigned int i,
  const double lambda
) {
  double cval = 0;
  CostResult cost_result;
  if (!contain(FASTCPD_FAMILIES, family)) {
    cost_result = get_optimized_cost(segment_start, segment_end);
  } else {
    if (warm_start && segment_end + 1 - segment_start >= 10 * p) {
      cost_result = get_cost_result(
        segment_start, segment_end, R_NilValue, lambda, false,
        wrap(
          segment_theta_hat.row(index_max(find(segment_indices <= segment_end))
        ).t())
        // Or use `wrap(start.col(segment_start))` for warm start.
      );
      update_start(segment_start, colvec(cost_result.par));
    } else {
      cost_result = get_cost_result(
        segment_start, segment_end, R_NilValue, lambda, false, R_NilValue
      );
    }
  }
  cval = cost_result.value;

  // If `vanilla_percentage` is not 1, then we need to keep track of
  // thetas for later `fastcpd` steps.
  if (
    vanilla_percentage < 1 && segment_end < vanilla_percentage * data_n_rows
  ) {
    update_theta_hat(i, cost_result.par);
    update_theta_sum(i, cost_result.par);
  }
  return cval;
}

double Fastcpd::get_cval_sen(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const unsigned int i,
  const double lambda
) {
  const unsigned int segment_length = segment_end - segment_start + 1;
  double cval = 0;
  update_cost_parameters(
    segment_end + 1, segment_start, i, k.get(), lambda, line_search
  );
  colvec theta = theta_sum.col(i) / segment_length;
  DEBUG_RCOUT(theta);
  if (!contain(FASTCPD_FAMILIES, family)) {
    Function cost_non_null = cost.get();
    SEXP cost_result = cost_non_null(
      data.rows(segment_start, segment_end), theta
    );
    cval = as<double>(cost_result);
  } else if (
    (family != "lasso" && segment_length >= p) ||
    (family == "lasso" && segment_length >= 3)
  ) {
    cval = get_cost_result(
      segment_start, segment_end, wrap(theta), lambda, false, R_NilValue
    ).value;
  }
  // else segment_length < p or for lasso segment_length < 3
  return cval;
}

CostResult Fastcpd::get_optimized_cost(
  const unsigned int segment_start,
  const unsigned int segment_end
) {
  Function cost_ = cost.get();
  CostResult cost_result;
  const mat data_segment = data.rows(segment_start, segment_end);
  if (cost_gradient.isNull() && cost_hessian.isNull()) {
    cost_result = {{colvec()}, {colvec()}, as<double>(cost_(data_segment))};
  } else if (p == 1) {
    Environment stats = Environment::namespace_env("stats");
    Function optim = stats["optim"];
    List optim_result = optim(
      Named("par") = 0,
      Named("fn") = InternalFunction(
        +[](double theta, mat data, Function cost_) {
          return cost_(
            Named("data") = data,
            Named("theta") = std::log(theta / (1 - theta))
          );
        }
      ),
      Named("method") = "Brent",
      Named("lower") = 0,
      Named("upper") = 1,
      Named("data") = data_segment,
      Named("cost") = cost_
    );
    colvec par = as<colvec>(optim_result["par"]);
    double value = as<double>(optim_result["value"]);
    cost_result =
      {{log(par / (1 - par))}, {colvec()}, exp(value) / (1 + exp(value))};
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
    cost_result =
      {{as<colvec>(optim_result["par"])}, {colvec()}, optim_result["value"]};
  } else {
    // # nocov start
    stop("This branch should not be reached at classes.cc: 945.");
    // # nocov end
  }
  return cost_result;
}

void Fastcpd::update_cost_parameters(
  const unsigned int t,
  const unsigned int tau,
  const unsigned int i,
  Function k,
  const double lambda,
  const colvec& line_search
) {
  List cost_update_result = update_cost_parameters_steps(
    0, t - 1, tau, i, k, momentum, lambda, line_search
  );
  DEBUG_RCOUT(line_search);
  update_theta_hat(i, as<colvec>(cost_update_result[0]));
  create_theta_sum(i, as<colvec>(cost_update_result[1]));
  update_hessian(i, as<mat>(cost_update_result[2]));
  update_momentum(as<colvec>(cost_update_result[3]));
}

void Fastcpd::update_cost_parameters_step(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const int i,
  const int data_start,
  const int data_end,
  const double lambda,
  const colvec& line_search
) {
  mat hessian_i = hessian.slice(i);
  colvec gradient;

  if (!contain(FASTCPD_FAMILIES, family)) {
    mat cost_hessian_result = cost_hessian_wrapper(
      segment_start + data_start, segment_start + data_end, theta_hat.col(i)
    );
    hessian_i += cost_hessian_result;
    colvec cost_gradient_result = cost_gradient_wrapper(
      segment_start + data_start, segment_start + data_end, theta_hat.col(i)
    );
    gradient = cost_gradient_result;
  } else {
    hessian_i += (this->*get_hessian)(
      segment_start + data_start, segment_start + data_end, theta_hat.col(i)
    );
    gradient = (this->*get_gradient)(
      segment_start + data_start, segment_start + data_end, theta_hat.col(i)
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
      colvec theta_candidate =
        theta_hat.col(i) + line_search[line_search_index] * momentum;
      colvec theta_upper_bound = arma::min(std::move(theta_candidate), upper);
      colvec theta_projected = arma::max(std::move(theta_upper_bound), lower);
      line_search_costs[line_search_index] = cost_function_wrapper(
        segment_start,
        segment_end,
        wrap(theta_projected),
        lambda,
        false,
        R_NilValue
      ).value;
    }
  }
  best_learning_rate = line_search[line_search_costs.index_min()];

  // Update theta_hat with momentum
  theta_hat.col(i) += best_learning_rate * momentum;

  theta_hat.col(i) = arma::min(theta_hat.col(i), upper);
  theta_hat.col(i) = arma::max(theta_hat.col(i), lower);

  if (family == "lasso" || family == "gaussian") {
    // Update theta_hat with L1 penalty
    double hessian_norm = norm(hessian_i, "fro");
    vec normd = abs(theta_hat.col(i)) - lambda / hessian_norm;
    theta_hat.col(i) = sign(theta_hat.col(i)) % arma::max(
      normd, zeros<colvec>(normd.n_elem)
    );
  }

  hessian.slice(i) = std::move(hessian_i);
}

List Fastcpd::update_cost_parameters_steps(
    const unsigned int segment_start,
    const unsigned int segment_end,
    const int tau,
    const int i,
    Function k,
    colvec momentum,
    const double lambda,
    const colvec& line_search
) {
  update_cost_parameters_step(
    segment_start,
    segment_end,
    i,
    0,
    segment_end - segment_start,
    lambda,
    line_search
  );

  const unsigned int segment_length = segment_end - segment_start + 1;

  for (int kk = 1; kk <= as<int>(k(segment_length - tau)); kk++) {
    for (unsigned j = tau + 1; j <= segment_length; j++) {
      update_cost_parameters_step(
        segment_start, segment_end, i, tau, j - 1, lambda, line_search
      );
    }
  }

  theta_sum.col(i) += theta_hat.col(i);
  return List::create(
    theta_hat.col(i), theta_sum.col(i), hessian.slice(i), momentum
  );
}

double Fastcpd::update_cost_value(
  double value,
  const unsigned int nrows
) {
  if (cost_adjustment == "MDL") {
    value = value * std::log2(M_E);
  }
  return value + get_cost_adjustment_value(nrows);
}

colvec Fastcpd::update_cp_set(const colvec raw_cp_set) {
  // Remove change points close to the boundaries.
  colvec cp_set = raw_cp_set;
  cp_set = cp_set(find(cp_set > trim * data_n_rows));
  cp_set = cp_set(find(cp_set < (1 - trim) * data_n_rows));
  colvec cp_set_ = zeros<vec>(cp_set.n_elem + 1);
  if (cp_set.n_elem) {
    cp_set_.rows(1, cp_set_.n_elem - 1) = std::move(cp_set);
  }
  cp_set = sort(unique(std::move(cp_set_)));

  // Remove change points close to each other.
  ucolvec cp_set_too_close = find(diff(cp_set) <= trim * data_n_rows);
  if (cp_set_too_close.n_elem > 0) {
    int rest_element_count = cp_set.n_elem - cp_set_too_close.n_elem;
    colvec cp_set_rest_left = zeros<vec>(rest_element_count),
          cp_set_rest_right = zeros<vec>(rest_element_count);
    for (unsigned int i = 0, i_left = 0, i_right = 0; i < cp_set.n_elem; i++) {
      if (
        ucolvec left_find = find(cp_set_too_close == i);
        left_find.n_elem == 0
      ) {
        cp_set_rest_left(i_left) = cp_set(i);
        i_left++;
      }
      if (
        ucolvec right_find = find(cp_set_too_close == i - 1);
        right_find.n_elem == 0
      ) {
        cp_set_rest_right(i_right) = cp_set(i);
        i_right++;
      }
    }
    cp_set = floor((cp_set_rest_left + cp_set_rest_right) / 2);
  }
  return cp_set(find(cp_set > 0));
}

void Fastcpd::update_err_sd(
  const unsigned int segment_index, const double err_var
) {
  err_sd(segment_index) = sqrt(err_var);
}

void Fastcpd::update_fastcpd_parameters(const unsigned int t) {
  const int segment_index = index_max(find(segment_indices <= t - 1));
  rowvec cum_coef_add = segment_theta_hat.row(segment_index),
             coef_add = segment_theta_hat.row(segment_index);
  mat hessian_new;
  if (family == "binomial") {
    const rowvec x = data.row(t - 1).tail(p);
    const double prob = 1 / (1 + exp(-dot(coef_add, x)));
    hessian_new = (x.t() * x) * prob * (1 - prob);
  } else if (family == "poisson") {
    const rowvec x = data.row(t - 1).tail(p);
    cum_coef_add = coef_add;
    hessian_new = (x.t() * x) * exp(dot(coef_add, x));
  } else if (family == "lasso" || family == "gaussian") {
    const rowvec x = data.row(t - 1).tail(p);
    hessian_new = x.t() * x + epsilon * eye<mat>(p, p);
  } else if (family == "arma") {
    hessian_new = (this->*get_hessian)(0, t - 1, coef_add.t());
  } else if (!contain(FASTCPD_FAMILIES, family)) {
    hessian_new = zeros<mat>(p, p);
  }

  update_theta_hat(coef_add.t());
  update_theta_sum(cum_coef_add.t());
  update_hessian(hessian_new);
}

void Fastcpd::update_hessian(
  const unsigned int slice, mat new_hessian
) {
  hessian.slice(slice) = new_hessian;
}

void Fastcpd::update_hessian(mat new_hessian) {
  hessian = join_slices(hessian, new_hessian);
}

void Fastcpd::update_hessian(ucolvec pruned_left) {
  hessian = hessian.slices(pruned_left);
}

void Fastcpd::update_momentum(colvec new_momentum) {
  momentum = new_momentum;
}

void Fastcpd::update_r_clock_tick(const std::string name) {
  if (!r_clock.empty()) {
    rClock.tick(name);
  }
}

void Fastcpd::update_r_clock_tock(const std::string name) {
  if (!r_clock.empty()) {
    rClock.tock(name);
  }
}

void Fastcpd::update_r_progress_start() {
  if (r_progress) {
    rProgress->tick(0);
  }
}

void Fastcpd::update_r_progress_tick() {
  if (r_progress) {
    rProgress->tick();
  }
}

void Fastcpd::update_start(const unsigned int col, const colvec start_col) {
  start.col(col) = start_col;
}

void Fastcpd::update_theta_hat(colvec new_theta_hat) {
  theta_hat = join_rows(theta_hat, new_theta_hat);
}

void Fastcpd::update_theta_hat(
  const unsigned int col, colvec new_theta_hat
) {
  theta_hat.col(col) = new_theta_hat;
}

void Fastcpd::update_theta_hat(ucolvec pruned_left) {
  theta_hat = theta_hat.cols(pruned_left);
}

void Fastcpd::update_theta_sum(colvec new_theta_sum) {
  theta_sum = join_rows(theta_sum, new_theta_sum);
}

void Fastcpd::update_theta_sum(ucolvec pruned_left) {
  theta_sum = theta_sum.cols(pruned_left);
}

}  // namespace fastcpd::classes
