#include "fastcpd_classes.h"
#include "fastcpd_constants.h"
#include "fastcpd_functions.h"
#include "RProgress.h"

using ::fastcpd::functions::negative_log_likelihood_arma;
using ::fastcpd::functions::negative_log_likelihood_glm;
using ::fastcpd::functions::negative_log_likelihood_lasso_cv;
using ::fastcpd::functions::negative_log_likelihood_lasso_wo_cv;
using ::fastcpd::functions::negative_log_likelihood_mean;
using ::fastcpd::functions::negative_log_likelihood_meanvariance;
using ::fastcpd::functions::negative_log_likelihood_mgaussian;
using ::fastcpd::functions::negative_log_likelihood_variance;

namespace fastcpd::classes {

Fastcpd::Fastcpd(
    const double beta,
    const double convexity_coef,
    Nullable<Function> cost,
    const string cost_adjustment,
    Nullable<Function> cost_gradient,
    Nullable<Function> cost_hessian,
    const bool cp_only,
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
    const bool r_progress,
    const int segment_count,
    const double trim,
    const colvec upper,
    const double vanilla_percentage,
    const mat variance_estimate,
    const bool warm_start
) : beta(beta),
    convexity_coef(convexity_coef),
    cost(cost),
    cost_adjustment(cost_adjustment),
    cost_gradient(cost_gradient),
    cost_hessian(cost_hessian),
    cp_only(cp_only),
    data(data),
    epsilon(epsilon),
    family(family),
    k(k),
    line_search(line_search),
    lower(lower),
    momentum_coef(momentum_coef),
    order(order),
    p(p),
    p_response(p_response),
    r_progress(r_progress),
    segment_count(segment_count),
    trim(trim),
    upper(upper),
    vanilla_percentage(vanilla_percentage),
    variance_estimate(variance_estimate),
    warm_start(warm_start) {
  n = data.n_rows;
  segment_indices = vec(n);
  segment_theta_hat = mat(segment_count, p);
  err_sd = vec(segment_count);
  act_num = vec(segment_count);
  start = zeros<mat>(p, n);
  theta_hat = mat(p, 1);
  theta_sum = mat(p, 1);
  hessian = cube(p, p, 1);
  momentum = vec(p);
  variance_data_mean = mean(data, 0);

  create_cost_function_wrapper(cost);
  create_cost_gradient_wrapper(cost_gradient);
  create_cost_hessian_wrapper(cost_hessian);

  // TODO(doccstat): Store environment functions from R.
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
    mat reversed_data = reverse(data, 0);
    colvec reversed_theta = reverse(theta);
    if (data.n_rows < max(order) + 1) {
      return ones(theta.n_elem);
    }
    colvec variance_term = zeros(data.n_rows);
    for (unsigned int i = max(order); i < data.n_rows; i++) {
      variance_term(i) = data(i, 0) - dot(
          reversed_theta.rows(order(1) + 1, sum(order)),
          data.rows(i - order(0), i - 1)
        ) - dot(
          reversed_theta.rows(1, order(1)),
          variance_term.rows(i - order(1), i - 1)
        );
    }
    colvec reversed_variance_term = reverse(variance_term);
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
    // Prevent numerical issues if `prob` is too large.
    hessian = (x.t() * x) * std::min(as_scalar(prob), 1e10);
  } else if (family == "lasso" || family == "gaussian") {
    hessian = x.t() * x;
  } else if (family == "arma") {
    // TODO(doccstat): Maybe we can store all these computations
    mat reversed_data = reverse(data, 0);
    colvec reversed_theta = reverse(theta);
    if (data.n_rows < max(order) + 1) {
      return eye(theta.n_elem, theta.n_elem);
    }
    colvec variance_term = zeros(data.n_rows);
    for (unsigned int i = max(order); i < data.n_rows; i++) {
      variance_term(i) = data(i, 0) - dot(
          reversed_theta.rows(order(1) + 1, sum(order)),
          data.rows(i - order(0), i - 1)
        ) - dot(
          reversed_theta.rows(1, order(1)),
          variance_term.rows(i - order(1), i - 1)
        );
    }
    colvec reversed_variance_term = reverse(variance_term);
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
    mat reversed_coef_phi = reverse(phi_coefficient, 0),
        reversed_coef_psi = reverse(psi_coefficient, 0);
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

void Fastcpd::create_theta_sum(
  const unsigned int col, colvec new_theta_sum
) {
  theta_sum.col(col) = new_theta_sum;
}

mat Fastcpd::get_theta_sum() {
  return theta_sum;
}

CostResult Fastcpd::negative_log_likelihood_wo_theta(
    mat data,
    double lambda,
    bool cv,
    Nullable<colvec> start
) {
  CostResult cost_result;
  if (family == "lasso" && cv) {
    cost_result = negative_log_likelihood_lasso_cv(data);
  } else if (family == "lasso" && !cv) {
    cost_result = negative_log_likelihood_lasso_wo_cv(data, lambda);
  } else if (
    family == "binomial" || family == "poisson" || family == "gaussian"
  ) {
    cost_result = negative_log_likelihood_glm(data, start, family);
  } else if (family == "arma") {
    cost_result = negative_log_likelihood_arma(data, order);
  } else if (family == "mean") {
    cost_result = negative_log_likelihood_mean(data, variance_estimate);
  } else if (family == "variance") {
    cost_result = negative_log_likelihood_variance(data, variance_data_mean);
  } else if (family == "meanvariance" || family == "mv") {
    cost_result = negative_log_likelihood_meanvariance(data, epsilon);
  } else if (family == "mgaussian") {
    cost_result = negative_log_likelihood_mgaussian(
      data, p_response, variance_estimate
    );
  } else {
    // # nocov start
    stop("This branch should not be reached at fastcpd_class_cost.cc: 193.");
    // # nocov end
  }
  return cost_result;
}

double Fastcpd::negative_log_likelihood_wo_cv(
    mat data,
    colvec theta,
    double lambda
) {
  vec y = data.col(0);
  if (family == "lasso" || family == "gaussian") {
    // Calculate negative log likelihood in gaussian family
    double penalty = lambda * accu(abs(theta));
    mat x = data.cols(1, data.n_cols - 1);
    return accu(square(y - x * theta)) / 2 + penalty;
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
    colvec reversed_theta = reverse(theta);
    if (data.n_rows < max(order) + 1) {
      return 0;
    }
    colvec variance_term = zeros(data.n_rows);
    for (unsigned int i = max(order); i < data.n_rows; i++) {
      variance_term(i) = data(i, 0) - dot(
          reversed_theta.rows(order(1) + 1, sum(order)),
          data.rows(i - order(0), i - 1)
        ) - dot(
          reversed_theta.rows(1, order(1)),
          variance_term.rows(i - order(1), i - 1)
        );
    }
    return (std::log(2.0 * M_PI) +
      std::log(theta(sum(order)))) * (data.n_rows - 2) / 2.0 +
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

  ucolvec r_t_set = {0};
  DEBUG_RCOUT(r_t_set);

  std::vector<colvec> cp_sets = {zeros<vec>(0)};
  linspace(1, n, n).for_each([&](int i) {
    cp_sets.push_back(zeros<vec>(1));
  });

  colvec fvec = zeros<vec>(n + 1);
  fvec.fill(arma::datum::inf);
  fvec(0) = -beta;
  DEBUG_RCOUT(fvec(0));

  RProgress::RProgress rProgress("[:bar] :current/:total in :elapsed", n);

  if (r_progress) {
    rProgress.tick(0);
  }

  if (contain(FASTCPD_FAMILIES, family) || vanilla_percentage < 1) {
    create_segment_indices();
    create_segment_statistics();
  }

  DEBUG_RCOUT(n);

  if (vanilla_percentage < 1) {
    create_gradients();
  }

  checkUserInterrupt();
  if (r_progress) {
    rProgress.tick();
  }

  for (int t = 1; t <= n; t++) {
    DEBUG_RCOUT(t);
    unsigned int r_t_count = r_t_set.n_elem;
    DEBUG_RCOUT(r_t_count);

    // Number of cost values is the same as the number of elements in R_t.
    colvec cval = zeros<vec>(r_t_count);

    // For tau in R_t \ {t-1}.
    for (unsigned int i = 1; i < r_t_count; i++) {
      cval(i - 1) = get_cval_for_r_t_set(r_t_set, i, t, lambda);
    }

    DEBUG_RCOUT(cval);
    cval(r_t_count - 1) = 0;

    if (vanilla_percentage != 1) {
      update_fastcpd_parameters(t);
    }

    // `beta` adjustment seems to work but there might be better choices.
    colvec obj = cval + fvec.rows(r_t_set) + beta;
    double min_obj = min(obj);
    double tau_star = r_t_set(index_min(obj));

    // Step 4
    cp_sets[t] = join_cols(cp_sets[tau_star], colvec{tau_star});
    DEBUG_RCOUT(cp_sets[t]);

    // Pruning step.
    ucolvec pruned_left =
      find(cval + fvec.rows(r_t_set) + convexity_coef <= min_obj);
    DEBUG_RCOUT(pruned_left);
    ucolvec pruned_r_t_set = zeros<ucolvec>(pruned_left.n_elem + 1);
    DEBUG_RCOUT(pruned_r_t_set);
    if (pruned_left.n_elem) {
      pruned_r_t_set.rows(0, pruned_left.n_elem - 1) = r_t_set(pruned_left);
    }
    DEBUG_RCOUT(pruned_r_t_set);
    pruned_r_t_set(pruned_left.n_elem) = t;
    r_t_set = std::move(pruned_r_t_set);
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
    if (r_progress) {
      rProgress.tick();
    }
  }

  return get_cp_set(cp_sets[n], lambda);
}

void Fastcpd::update_theta_sum(
  const unsigned int col, colvec new_theta_sum
) {
  theta_sum.col(col) += new_theta_sum;
}

}  // namespace fastcpd::classes
