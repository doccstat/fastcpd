#include "fastcpd_class.h"
#include "RProgress.h"

using ::arma::abs;
using ::arma::accu;
using ::arma::as_scalar;
using ::arma::diff;
using ::arma::dot;
using ::arma::find;
using ::arma::eye;
using ::arma::floor;
using ::arma::index_max;
using ::arma::index_min;
using ::arma::join_cols;
using ::arma::join_rows;
using ::arma::join_slices;
using ::arma::linspace;
using ::arma::max;
using ::arma::mean;
using ::arma::min;
using ::arma::norm;
using ::arma::sign;
using ::arma::solve;
using ::arma::sort;
using ::arma::square;
using ::arma::unique;
using ::arma::zeros;
using ::Rcpp::as;
using ::Rcpp::checkUserInterrupt;
using ::Rcpp::Named;
using ::Rcpp::wrap;
using ::std::make_unique;

using ::arma::vec;
using ::Rcpp::Environment;
using ::Rcpp::InternalFunction;

using ::Rcpp::Rcout;

#define ERROR(msg) \
  Rcout << "error: " << __FILE__ << ": " << \
                        __LINE__ << ": " << msg << std::endl
#define FATAL(msg) \
  Rcout << "fatal: " << __FILE__ << ": " << \
                        __LINE__ << ": " << msg << std::endl; \
  throw std::runtime_error(msg)

namespace fastcpd::classes {

Fastcpd::Fastcpd(
    const double beta,
    Nullable<Function> cost,
    const string cost_adjustment,
    Nullable<Function> cost_gradient,
    Nullable<Function> cost_hessian,
    const bool cp_only,
    const unsigned int d,
    mat data,
    const double epsilon,
    const string family,
    Nullable<Function> k,
    colvec line_search,
    const colvec lower,
    const double momentum_coef,
    const colvec order,
    const int p,
    const unsigned int p_response,
    const double pruning_coef,
    const string r_clock,
    const bool r_progress,
    const int segment_count,
    const double trim,
    const colvec upper,
    const double vanilla_percentage,
    const mat variance_estimate,
    const bool warm_start
) : act_num(colvec(segment_count)),
    beta(beta),
    cost_adjustment(cost_adjustment),
    cp_only(cp_only),
    d(d),
    data(data),
    data_n_rows(data.n_rows),
    data_n_cols(data.n_cols),
    epsilon(epsilon),
    err_sd(colvec(segment_count)),
    family(family),
    hessian(cube(p, p, 1)),
    line_search(line_search),
    lower(lower),
    momentum(vec(p)),
    momentum_coef(momentum_coef),
    order(order),
    p(p),
    p_response(p_response),
    pruning_coef(pruning_coef),
    r_clock(r_clock),
    r_progress(r_progress),
    rProgress(make_unique<RProgress::RProgress>(kRProgress, data_n_rows)),
    segment_count(segment_count),
    segment_indices(round(linspace(0, data_n_rows, segment_count + 1))),
    segment_theta_hat(mat(segment_count, p)),
    start(zeros<mat>(p, data_n_rows)),
    theta_hat(mat(p, 1)),
    theta_sum(mat(p, 1)),
    trim(trim),
    upper(upper),
    vanilla_percentage(vanilla_percentage),
    variance_estimate(variance_estimate),
    warm_start(warm_start),
    zero_data(join_cols(zeros<rowvec>(data_n_cols), data)) {
  if (k.isNotNull()) {
    this->k = std::make_unique<Function>(k);
  }

  zero_data_c = (double **)malloc((data_n_rows + 1) * sizeof(double *));
  if (zero_data_c == NULL) {
    FATAL("Memory allocation failed.");  // # nocov
  }
  unsigned int i = 0;
  while (i < data_n_rows + 1) {
    zero_data_c[i] = (double *)malloc(data_n_cols * sizeof(double));
    if (zero_data_c[i] == NULL) {
      break;  // # nocov
    }
    i++;
  }
  if (i < data_n_rows + 1) {
    for (unsigned int j = 0; j < i; j++) {  // # nocov start
      free(zero_data_c[j]);
    }
    free(zero_data_c);
    FATAL("Memory allocation failed.");  // # nocov end
  }
  for (unsigned int i = 0; i < data_n_rows + 1; i++) {
    for (unsigned int j = 0; j < data_n_cols; j++) {
      zero_data_c[i][j] = zero_data(i, j);
    }
  }

  create_gets(cost, cost_gradient, cost_hessian);
}

List Fastcpd::run() {
  // Set up the initial values.
  double lambda = 0;

  ucolvec r_t_set = zeros<ucolvec>(data_n_rows);
  r_t_set(1) = 1;
  unsigned int r_t_count = 2;

  std::vector<colvec> cp_sets = {{0}};
  linspace(1, data_n_rows, data_n_rows).for_each([&](int i) {
    cp_sets.push_back({0});
  });

  colvec fvec = zeros<vec>(data_n_rows + 1);
  fvec.fill(arma::datum::inf);
  fvec(0) = -beta;

  create_segment_statistics();
  create_gradients();

  checkUserInterrupt();
  update_r_progress_start();
  update_r_progress_tick();

  for (unsigned int t = 2; t <= data_n_rows; t++) {
    colvec cval = zeros<vec>(r_t_count);

    update_r_clock_tick("r_t_set_for_loop");
    for (unsigned int i = 0; i < r_t_count - 1; i++) {
      cval(i) = get_cval_for_r_t_set(r_t_set(i), i, t, lambda);
    }
    update_r_clock_tock("r_t_set_for_loop");
    update_r_clock_tick("pruning");

    if (vanilla_percentage != 1) {
      update_fastcpd_parameters(t);
    }

    colvec obj = cval + fvec.rows(r_t_set.rows(0, r_t_count - 1)) + beta;
    double min_obj = min(obj);
    double tau_star = r_t_set(index_min(obj));

    cp_sets[t] = join_cols(cp_sets[tau_star], colvec{tau_star});

    ucolvec pruned_left = find(
      cval + fvec.rows(r_t_set.rows(0, r_t_count - 1)) + pruning_coef <= min_obj
    );
    r_t_count = pruned_left.n_elem + 1;
    if (pruned_left.n_elem) {
      r_t_set.rows(0, pruned_left.n_elem - 1) = r_t_set(pruned_left);
    }
    r_t_set(pruned_left.n_elem) = t;

    if (vanilla_percentage != 1) {
      update_theta_hat(pruned_left);
      update_theta_sum(pruned_left);
      update_hessian(pruned_left);
    }

    fvec(t) = min_obj;

    checkUserInterrupt();
    update_r_progress_tick();
    update_r_clock_tock("pruning");
  }

  create_clock_in_r(r_clock);

  List result = get_cp_set(cp_sets[data_n_rows], lambda);

  for (unsigned int i = 0; i < data_n_rows + 1; i++) {
    free(zero_data_c[i]);
  }
  free(zero_data_c);

  return result;
}

void Fastcpd::create_clock_in_r(const std::string name) {
  if (!r_clock.empty()) {
    rClock.stop(name);
  }
}

void Fastcpd::create_gets(
  Nullable<Function>& cost,
  Nullable<Function>& cost_gradient,
  Nullable<Function>& cost_hessian
) {
  if (family == "arma" && order(0) > 0) {
    get_gradient = &Fastcpd::get_gradient_arma;
    get_hessian = &Fastcpd::get_hessian_arma;
    get_nll_sen = &Fastcpd::get_nll_sen_arma;
    get_nll_pelt = &Fastcpd::get_nll_pelt_arma;
  } else if (family.compare("binomial") == 0) {
    get_gradient = &Fastcpd::get_gradient_binomial;
    get_hessian = &Fastcpd::get_hessian_binomial;
    get_nll_sen = &Fastcpd::get_nll_sen_binomial;
    get_nll_pelt = &Fastcpd::get_nll_pelt_glm;
  } else if (family == "gaussian") {
    get_gradient = &Fastcpd::get_gradient_lm;
    get_hessian = &Fastcpd::get_hessian_lm;
    get_nll_sen = &Fastcpd::get_nll_sen_lm;
    get_nll_pelt = &Fastcpd::get_nll_pelt_glm;
  } else if (family == "lasso") {
    get_gradient = &Fastcpd::get_gradient_lm;
    get_hessian = &Fastcpd::get_hessian_lm;
    get_nll_sen = &Fastcpd::get_nll_sen_lm;
    get_nll_pelt = &Fastcpd::get_nll_pelt_lasso;
  } else if (family == "arma" && order(0) == 0) {
    get_gradient = &Fastcpd::get_gradient_ma;
    get_hessian = &Fastcpd::get_hessian_ma;
    get_nll_sen = &Fastcpd::get_nll_sen_ma;
    get_nll_pelt = &Fastcpd::get_nll_pelt_arma;
  } else if (family == "mean") {
    get_nll_pelt = &Fastcpd::get_nll_pelt_mean;
  } else if (family == "meanvariance") {
    get_nll_pelt = &Fastcpd::get_nll_pelt_meanvariance;
  } else if (family == "mgaussian") {
    get_nll_pelt = &Fastcpd::get_nll_pelt_mgaussian;
  } else if (family.compare("poisson") == 0) {
    get_gradient = &Fastcpd::get_gradient_poisson;
    get_hessian = &Fastcpd::get_hessian_poisson;
    get_nll_sen = &Fastcpd::get_nll_sen_poisson;
    get_nll_pelt = &Fastcpd::get_nll_pelt_glm;
  } else if (family == "variance") {
    get_nll_pelt = &Fastcpd::get_nll_pelt_variance;
  } else {
    this->cost = make_unique<Function>(cost);
    if (cost_gradient.isNotNull() || cost_hessian.isNotNull()) {
      this->cost_gradient = make_unique<Function>(cost_gradient);
      this->cost_hessian = make_unique<Function>(cost_hessian);
    }
    get_gradient = &Fastcpd::get_gradient_custom;
    get_hessian = &Fastcpd::get_hessian_custom;
    get_nll_sen = &Fastcpd::get_nll_sen_custom;
    get_nll_pelt = &Fastcpd::get_nll_pelt_custom;
  }

  // TODO(doccstat): Store environment functions from R.
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
  } else if (family == "custom") {
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    theta_sum.col(0) = segment_theta_hat.row(0).t();
    hessian.slice(0) = zeros<mat>(p, p);
  }
}

// TODO(doccstat): Use `segment_theta` as warm start.

void Fastcpd::create_segment_statistics() {
  if (family == "custom" && vanilla_percentage == 1) return;
  for (
    int segment_index = 0; segment_index < segment_count; ++segment_index
  ) {
    rowvec segment_theta = get_cost_result(
      segment_indices(segment_index),
      segment_indices(segment_index + 1) - 1,
      R_NilValue,
      0,
      true,
      R_NilValue
    ).par;

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

void Fastcpd::create_theta_sum(
  const unsigned int col, colvec new_theta_sum
) {
  theta_sum.col(col) = new_theta_sum;
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
    cost_result = (this->*get_nll_pelt)(
      segment_start, segment_end, lambda, cv, start
    );
  } else {
    cost_result = CostResult{
      {colvec()},
      {colvec()},
      (this->*get_nll_sen)(
        segment_start, segment_end, as<colvec>(theta), lambda
      )
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
    CostResult cost_result = get_cost_result(
      cp_loc(i), cp_loc(i + 1) - 1, R_NilValue, lambda, false, R_NilValue
    );

    cost_values(i) = cost_result.value;

    // Parameters are not involved for PELT.
    if (vanilla_percentage < 1) {
      thetas.col(i) = colvec(cost_result.par);
    }

    // Residual is only calculated for built-in families.
    if (
      family != "custom" && family != "mean" &&
      family != "variance" && family != "meanvariance"
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
  if (
    (family == "binomial" || family == "poisson") &&
    (warm_start && segment_end + 1 - segment_start >= 10 * p)
  ) {
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
    segment_end + 1, segment_start, i, lambda, line_search
  );
  colvec theta = theta_sum.col(i) / segment_length;
  if (family == "custom") {
    cval = (this->*get_nll_sen)(
      segment_start, segment_end, theta, lambda
    );
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
  CostResult cost_result;
  const mat data_segment = data.rows(segment_start, segment_end);
  if (p == 1) {
    Environment stats = Environment::namespace_env("stats");
    Function optim = stats["optim"];
    List optim_result = optim(
      Named("par") = 0,
      Named("fn") = InternalFunction(
        +[](double theta, mat data, Function cost) {
          return cost(
            Named("data") = data,
            Named("theta") = std::log(theta / (1 - theta))
          );
        }
      ),
      Named("method") = "Brent",
      Named("lower") = 0,
      Named("upper") = 1,
      Named("data") = data_segment,
      Named("cost") = *cost
    );
    colvec par = as<colvec>(optim_result["par"]);
    double value = as<double>(optim_result["value"]);
    cost_result =
      {{log(par / (1 - par))}, {colvec()}, exp(value) / (1 + exp(value))};
  } else {
    Environment stats = Environment::namespace_env("stats");
    Function optim = stats["optim"];
    List optim_result = optim(
      Named("par") = zeros<vec>(p),
      Named("fn") = *cost,
      Named("method") = "L-BFGS-B",
      Named("data") = data_segment,
      Named("lower") = lower,
      Named("upper") = upper
    );
    cost_result =
      {{as<colvec>(optim_result["par"])}, {colvec()}, optim_result["value"]};
  }
  return cost_result;
}

void Fastcpd::update_cost_parameters(
  const unsigned int t,
  const unsigned int tau,
  const unsigned int i,
  const double lambda,
  const colvec& line_search
) {
  List cost_update_result = update_cost_parameters_steps(
    0, t - 1, tau, i, momentum, lambda, line_search
  );
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

  hessian_i += (this->*get_hessian)(
    segment_start + data_start, segment_start + data_end, theta_hat.col(i)
  );
  gradient = (this->*get_gradient)(
    segment_start + data_start, segment_start + data_end, theta_hat.col(i)
  );

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
      line_search_costs[line_search_index] = (this->*get_nll_sen)(
        segment_start, segment_end, theta_projected, lambda
      );
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

  for (int kk = 1; kk <= as<int>((*k)(segment_length - tau)); kk++) {
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
  } else if (family == "custom") {
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

void Fastcpd::update_theta_sum(
  const unsigned int col, colvec new_theta_sum
) {
  theta_sum.col(col) += new_theta_sum;
}

}  // namespace fastcpd::classes
