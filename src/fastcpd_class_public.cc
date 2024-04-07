#include "fastcpd_classes.h"

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
    cost(cost),
    cost_adjustment(cost_adjustment),
    cost_gradient(cost_gradient),
    cost_hessian(cost_hessian),
    cp_only(cp_only),
    d(d),
    data(data),
    data_n_rows(data.n_rows),
    data_n_cols(data.n_cols),
    epsilon(epsilon),
    err_sd(colvec(segment_count)),
    family(family),
    hessian(cube(p, p, 1)),
    k(k),
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
  rProgress = std::make_unique<RProgress::RProgress>(
    "[:bar] :current/:total in :elapsed", data_n_rows
  );

  create_cost_function_wrapper(cost);
  create_cost_gradient_wrapper(cost_gradient);
  create_cost_hessian_wrapper(cost_hessian);

  // TODO(doccstat): Store environment functions from R.
}

colvec Fastcpd::cost_update_gradient(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  const unsigned int segment_length = segment_end - segment_start + 1;
  const mat data_segment = data.rows(segment_start, segment_end);
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  double y = new_data(0);
  colvec gradient;
  if (family.compare("binomial") == 0) {
    gradient = - (y - 1 / (1 + exp(-as_scalar(x * theta)))) * x.t();
  } else if (family.compare("poisson") == 0) {
    gradient = - (y - exp(as_scalar(x * theta))) * x.t();
  } else if (family == "lasso" || family == "gaussian") {
    gradient = - (y - as_scalar(x * theta)) * x.t();
  } else if (family == "arma" && order(0) == 0) {
    const unsigned int q = order(1);
    mat reversed_data = reverse(data_segment, 0);
    colvec reversed_theta = reverse(theta);
    if (segment_length < q + 1) {
      return ones(theta.n_elem);
    }
    colvec variance_term = zeros(segment_length);
    for (unsigned int i = q; i < segment_length; i++) {
      variance_term(i) = data_segment(i, 0) -
        dot(reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1));
    }
    colvec reversed_variance_term = reverse(variance_term);
    mat psi_coefficient = zeros(segment_length, q);
    for (unsigned int i = q; i < segment_length; i++) {
      psi_coefficient.row(i) = -reversed_variance_term.rows(
          segment_length - i, segment_length - i + q - 1
        ).t() - reversed_theta.rows(1, q).t() *
        psi_coefficient.rows(i - q, i - 1);
    }
    gradient = zeros(q + 1);
    gradient.rows(0, q - 1) =
      psi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(q);
    gradient(q) = 1.0 / 2.0 / theta(q) -
      std::pow(variance_term(segment_length - 1), 2) / 2.0 /
      std::pow(theta(q), 2);
  } else if (family == "arma" && order(0) > 0) {
    mat reversed_data = reverse(data_segment, 0);
    colvec reversed_theta = reverse(theta);
    if (segment_length < max(order) + 1) {
      return ones(theta.n_elem);
    }
    colvec variance_term = zeros(segment_length);
    for (unsigned int i = max(order); i < segment_length; i++) {
      variance_term(i) = data_segment(i, 0) - dot(
          reversed_theta.rows(order(1) + 1, sum(order)),
          data_segment.rows(i - order(0), i - 1)
        ) - dot(
          reversed_theta.rows(1, order(1)),
          variance_term.rows(i - order(1), i - 1)
        );
    }
    colvec reversed_variance_term = reverse(variance_term);
    mat phi_coefficient = zeros(segment_length, order(0)),
        psi_coefficient = zeros(segment_length, order(1));
    for (unsigned int i = max(order); i < segment_length; i++) {
      phi_coefficient.row(i) = -reversed_data.rows(
        segment_length - i, segment_length - i + order(0) - 1
      ).t() - reversed_theta.rows(1, order(1)).t() *
      phi_coefficient.rows(i - order(1), i - 1);
    }
    for (unsigned int i = order(1); i < segment_length; i++) {
      psi_coefficient.row(i) = -reversed_variance_term.rows(
          segment_length - i, segment_length - i + order(1) - 1
        ).t() - reversed_theta.rows(1, order(1)).t() *
        psi_coefficient.rows(i - order(1), i - 1);
    }
    gradient = zeros(sum(order) + 1);
    gradient.rows(0, order(0) - 1) =
      phi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(sum(order));
    gradient.rows(order(0), sum(order) - 1) =
      psi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(sum(order));
    gradient(sum(order)) = 1.0 / 2.0 / theta(sum(order)) -
      std::pow(variance_term(segment_length - 1), 2) / 2.0 /
      std::pow(theta(sum(order)), 2);
  } else {
    // # nocov start
    stop("This branch should not be reached at functions.cc: 188.");
    // # nocov end
  }
  return gradient;
}

mat Fastcpd::cost_update_hessian(
  const unsigned int segment_start,
  const unsigned int segment_end,
  const colvec& theta
) {
  const mat data_segment = data.rows(segment_start, segment_end);
  const unsigned int segment_length = segment_end - segment_start + 1;
  rowvec new_data = data_segment.row(segment_length - 1);
  rowvec x = new_data.tail(new_data.n_elem - 1);
  mat hessian;
  if (family.compare("binomial") == 0) {
    double prob = 1 / (1 + exp(-as_scalar(x * theta)));
    hessian = (x.t() * x) * as_scalar((1 - prob) * prob);
  } else if (family.compare("poisson") == 0) {
    double prob = exp(as_scalar(x * theta));
    // Prevent numerical issues if `prob` is too large.
    hessian = (x.t() * x) * std::min(as_scalar(prob), 1e10);
  } else if (family == "lasso" || family == "gaussian") {
    hessian = x.t() * x;
  } else if (family == "arma" && order(0) == 0) {
    const unsigned int q = order(1);
    // TODO(doccstat): Maybe we can store all these computations
    mat reversed_data = reverse(data_segment, 0);
    colvec reversed_theta = reverse(theta);
    if (segment_length < q + 1) {
      return eye(theta.n_elem, theta.n_elem);
    }
    colvec variance_term = zeros(segment_length);
    for (unsigned int i = q; i < segment_length; i++) {
      variance_term(i) = data_segment(i, 0) - dot(
        reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1)
      );
    }
    colvec reversed_variance_term = reverse(variance_term);
    mat psi_coefficient = zeros(segment_length, q);
    for (unsigned int i = q; i < segment_length; i++) {
      psi_coefficient.row(i) = -reversed_variance_term.rows(
          segment_length - i, segment_length - i + q - 1
        ).t() - reversed_theta.rows(1, q).t() *
        psi_coefficient.rows(i - q, i - 1);
    }
    mat reversed_coef_psi = reverse(psi_coefficient, 0);
    cube psi_psi_coefficient = zeros(q, q, segment_length);
    for (unsigned int i = q; i < segment_length; i++) {
      mat psi_psi_coefficient_part = zeros(q, q);
      for (unsigned int j = 1; j <= q; j++) {
        psi_psi_coefficient_part +=
          psi_psi_coefficient.slice(i - j) * theta(j - 1);
      }
      psi_psi_coefficient.slice(i) = -reversed_coef_psi.rows(
          segment_length - i, segment_length - i + q - 1
        ) - reversed_coef_psi.rows(
          segment_length - i, segment_length - i + q - 1
        ).t() - psi_psi_coefficient_part;
    }
    hessian = zeros(q + 1, q + 1);
    hessian.submat(0, 0, q - 1, q - 1) = (
      psi_coefficient.row(segment_length - 1).t() *
      psi_coefficient.row(segment_length - 1) +
      psi_psi_coefficient.slice(segment_length - 1) *
      variance_term(segment_length - 1)
    ) / theta(q);
    hessian.submat(0, q, q - 1, q) =
      -psi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(q) / theta(q);
    hessian.submat(q, 0, q, q - 1) = hessian.submat(0, q, q - 1, q).t();
    hessian(q, q) =
      std::pow(variance_term(segment_length - 1), 2) /
      std::pow(theta(q), 3) -
      1.0 / 2.0 / std::pow(theta(q), 2);
  } else if (family == "arma" && order(0) > 0) {
    // TODO(doccstat): Maybe we can store all these computations
    mat reversed_data = reverse(data_segment, 0);
    colvec reversed_theta = reverse(theta);
    if (segment_length < max(order) + 1) {
      return eye(theta.n_elem, theta.n_elem);
    }
    colvec variance_term = zeros(segment_length);
    for (unsigned int i = max(order); i < segment_length; i++) {
      variance_term(i) = data_segment(i, 0) - dot(
          reversed_theta.rows(order(1) + 1, sum(order)),
          data_segment.rows(i - order(0), i - 1)
        ) - dot(
          reversed_theta.rows(1, order(1)),
          variance_term.rows(i - order(1), i - 1)
        );
    }
    colvec reversed_variance_term = reverse(variance_term);
    mat phi_coefficient = zeros(segment_length, order(0)),
        psi_coefficient = zeros(segment_length, order(1));
    for (unsigned int i = max(order); i < segment_length; i++) {
      phi_coefficient.row(i) = -reversed_data.rows(
        segment_length - i, segment_length - i + order(0) - 1
      ).t() - reversed_theta.rows(1, order(1)).t() *
      phi_coefficient.rows(i - order(1), i - 1);
    }
    for (unsigned int i = order(1); i < segment_length; i++) {
      psi_coefficient.row(i) = -reversed_variance_term.rows(
          segment_length - i, segment_length - i + order(1) - 1
        ).t() - reversed_theta.rows(1, order(1)).t() *
        psi_coefficient.rows(i - order(1), i - 1);
    }
    mat reversed_coef_phi = reverse(phi_coefficient, 0),
        reversed_coef_psi = reverse(psi_coefficient, 0);
    cube phi_psi_coefficient = zeros(order(1), order(0), segment_length),
         psi_psi_coefficient = zeros(order(1), order(1), segment_length);
    for (unsigned int i = order(1); i < segment_length; i++) {
      mat phi_psi_coefficient_part = zeros(order(1), order(0)),
          psi_psi_coefficient_part = zeros(order(1), order(1));
      for (unsigned int j = 1; j <= order(1); j++) {
        phi_psi_coefficient_part +=
          phi_psi_coefficient.slice(i - j) * theta(order(0) - 1 + j);
      }
      phi_psi_coefficient.slice(i) = -reversed_coef_phi.rows(
        segment_length - i, segment_length - i + order(1) - 1
      ) - phi_psi_coefficient_part;
      for (unsigned int j = 1; j <= order(1); j++) {
        psi_psi_coefficient_part +=
          psi_psi_coefficient.slice(i - j) * theta(order(0) - 1 + j);
      }
      psi_psi_coefficient.slice(i) = -reversed_coef_psi.rows(
          segment_length - i, segment_length - i + order(1) - 1
        ) - reversed_coef_psi.rows(
          segment_length - i, segment_length - i + order(1) - 1
        ).t() - psi_psi_coefficient_part;
    }
    hessian = zeros(sum(order) + 1, sum(order) + 1);
    hessian.submat(0, 0, order(0) - 1, order(0) - 1) =
      phi_coefficient.row(segment_length - 1).t() *
      phi_coefficient.row(segment_length - 1) / theta(sum(order));
    hessian.submat(0, order(0), order(0) - 1, sum(order) - 1) = (
      phi_psi_coefficient.slice(segment_length - 1).t() *
        variance_term(segment_length - 1) +
        phi_coefficient.row(segment_length - 1).t() *
        psi_coefficient.row(segment_length - 1)
    ) / theta(sum(order));
    hessian.submat(order(0), 0, sum(order) - 1, order(0) - 1) =
      hessian.submat(0, order(0), order(0) - 1, sum(order) - 1).t();
    hessian.submat(0, sum(order), order(0) - 1, sum(order)) =
      -phi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(sum(order)) / theta(sum(order));
    hessian.submat(sum(order), 0, sum(order), order(0) - 1) =
      hessian.submat(0, sum(order), order(0) - 1, sum(order)).t();
    hessian.submat(order(0), order(0), sum(order) - 1, sum(order) - 1) = (
      psi_coefficient.row(segment_length - 1).t() *
      psi_coefficient.row(segment_length - 1) +
      psi_psi_coefficient.slice(segment_length - 1) *
      variance_term(segment_length - 1)
    ) / theta(sum(order));
    hessian.submat(order(0), sum(order), sum(order) - 1, sum(order)) =
      -psi_coefficient.row(segment_length - 1).t() *
      variance_term(segment_length - 1) / theta(sum(order)) / theta(sum(order));
    hessian.submat(sum(order), order(0), sum(order), sum(order) - 1) =
      hessian.submat(order(0), sum(order), sum(order) - 1, sum(order)).t();
    hessian(sum(order), sum(order)) =
      std::pow(variance_term(segment_length - 1), 2) /
      std::pow(theta(sum(order)), 3) -
      1.0 / 2.0 / std::pow(theta(sum(order)), 2);
  }
  return hessian;
}

void Fastcpd::create_theta_sum(
  const unsigned int col, colvec new_theta_sum
) {
  theta_sum.col(col) = new_theta_sum;
}

mat Fastcpd::get_theta_sum() {
  return theta_sum;
}

CostResult Fastcpd::get_nll_wo_theta(
    const unsigned int segment_start,
    const unsigned int segment_end,
    double lambda,
    bool cv,
    Nullable<colvec> start
) {
  CostResult cost_result;
  if (family == "lasso" && cv) {
    cost_result = get_nll_lasso_cv(segment_start, segment_end);
  } else if (family == "lasso" && !cv) {
    cost_result = get_nll_lasso_wo_cv(segment_start, segment_end, lambda);
  } else if (
    family == "binomial" || family == "poisson" || family == "gaussian"
  ) {
    cost_result = get_nll_glm(segment_start, segment_end, start);
  } else if (family == "arma") {
    cost_result = get_nll_arma(segment_start, segment_end);
  } else if (family == "mean") {
    cost_result = get_nll_mean(segment_start, segment_end);
  } else if (family == "variance") {
    cost_result = get_nll_variance(segment_start, segment_end);
  } else if (family == "meanvariance") {
    cost_result = get_nll_meanvariance(segment_start, segment_end);
  } else if (family == "mgaussian") {
    cost_result = get_nll_mgaussian(segment_start, segment_end);
  } else {
    // # nocov start
    stop("This branch should not be reached at fastcpd_class_cost.cc: 193.");
    // # nocov end
  }
  return cost_result;
}

double Fastcpd::get_nll_wo_cv(
    const unsigned int segment_start,
    const unsigned int segment_end,
    colvec theta,
    double lambda
) {
  mat data_segment = data.rows(segment_start, segment_end);
  vec y = data_segment.col(0);
  if (family == "lasso" || family == "gaussian") {
    // Calculate negative log likelihood in gaussian family
    double penalty = lambda * accu(abs(theta));
    mat x = data_segment.cols(1, data_segment.n_cols - 1);
    return accu(square(y - x * theta)) / 2 + penalty;
  } else if (family == "binomial") {
    // Calculate negative log likelihood in binomial family
    mat x = data_segment.cols(1, data_segment.n_cols - 1);
    colvec u = x * theta;
    return accu(-y % u + arma::log(1 + exp(u)));
  } else if (family == "poisson") {
    mat x = data_segment.cols(1, data_segment.n_cols - 1);
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
  } else if (family == "arma" && order(0) == 0) {
    const unsigned int q = order(1);
    colvec reversed_theta = reverse(theta);
    if (data_segment.n_rows < q + 1) {
      return 0;
    }
    colvec variance_term = zeros(data_segment.n_rows);
    for (unsigned int i = q; i < data_segment.n_rows; i++) {
      variance_term(i) = data_segment(i, 0) - dot(
          reversed_theta.rows(1, q), variance_term.rows(i - q, i - 1)
        );
    }
    return (std::log(2.0 * M_PI) +
      std::log(theta(q))) * (data_segment.n_rows - 2) / 2.0 +
      dot(variance_term, variance_term) / 2.0 / theta(q);
  } else if (family == "arma" && order(0) > 0) {
    colvec reversed_theta = reverse(theta);
    if (data_segment.n_rows < max(order) + 1) {
      return 0;
    }
    colvec variance_term = zeros(data_segment.n_rows);
    for (unsigned int i = max(order); i < data_segment.n_rows; i++) {
      variance_term(i) = data_segment(i, 0) - dot(
          reversed_theta.rows(order(1) + 1, sum(order)),
          data_segment.rows(i - order(0), i - 1)
        ) - dot(
          reversed_theta.rows(1, order(1)),
          variance_term.rows(i - order(1), i - 1)
        );
    }
    return (std::log(2.0 * M_PI) +
      std::log(theta(sum(order)))) * (data_segment.n_rows - 2) / 2.0 +
      dot(variance_term, variance_term) / 2.0 / theta(sum(order));
  } else {
    // # nocov start
    stop("This branch should not be reached at functions.cc: 248.");
    // # nocov end
  }
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
    DEBUG_RCOUT(t);
    DEBUG_RCOUT(r_t_count);

    colvec cval = zeros<vec>(r_t_count);

    update_r_clock_tick("r_t_set_for_loop");
    for (unsigned int i = 0; i < r_t_count - 1; i++) {
      cval(i) = get_cval_for_r_t_set(r_t_set(i), i, t, lambda);
    }
    update_r_clock_tock("r_t_set_for_loop");
    update_r_clock_tick("pruning");

    DEBUG_RCOUT(cval);

    if (vanilla_percentage != 1) {
      update_fastcpd_parameters(t);
    }

    colvec obj = cval + fvec.rows(r_t_set.rows(0, r_t_count - 1)) + beta;
    DEBUG_RCOUT(obj);
    double min_obj = min(obj);
    double tau_star = r_t_set(index_min(obj));
    DEBUG_RCOUT(tau_star);

    cp_sets[t] = join_cols(cp_sets[tau_star], colvec{tau_star});
    DEBUG_RCOUT(cp_sets[t]);

    ucolvec pruned_left = find(
      cval + fvec.rows(r_t_set.rows(0, r_t_count - 1)) + pruning_coef <= min_obj
    );
    r_t_count = pruned_left.n_elem + 1;
    if (pruned_left.n_elem) {
      r_t_set.rows(0, pruned_left.n_elem - 1) = r_t_set(pruned_left);
    }
    r_t_set(pruned_left.n_elem) = t;
    DEBUG_RCOUT(r_t_set);

    if (vanilla_percentage != 1) {
      update_theta_hat(pruned_left);
      update_theta_sum(pruned_left);
      update_hessian(pruned_left);
    }

    // Objective function F(t).
    fvec(t) = min_obj;
    DEBUG_RCOUT(fvec.rows(0, t));

    checkUserInterrupt();
    update_r_progress_tick();
    update_r_clock_tock("pruning");
  }

  create_clock_in_r(r_clock);

  return get_cp_set(cp_sets[data_n_rows], lambda);
}

void Fastcpd::update_theta_sum(
  const unsigned int col, colvec new_theta_sum
) {
  theta_sum.col(col) += new_theta_sum;
}

}  // namespace fastcpd::classes
