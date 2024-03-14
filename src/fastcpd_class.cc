#ifdef EXPERIMENT
#include <algorithm>
#include <execution>
#endif

#include "fastcpd_classes.h"
#include "fastcpd_constants.h"
#include "RProgress.h"

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
    const bool pruning,
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
    pruning(pruning),
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

double Fastcpd::adjust_cost_value(
  double value,
  const unsigned int nrows
) {
  if (cost_adjustment == "MDL") {
    value = value * std::log2(M_E);
  }
  return value + get_cost_adjustment_value(nrows);
}

double Fastcpd::get_cost_adjustment_value(const unsigned nrows) {
  double adjusted = 0;
  if (cost_adjustment == "MBIC" || cost_adjustment == "MDL") {
    adjusted = data.n_cols * std::log(nrows) / 2;
  }
  if (cost_adjustment == "MDL") {
    adjusted *= std::log2(M_E);
  }
  return adjusted;
}

List Fastcpd::get_optimized_cost(const mat data_segment) {
  Function cost_ = cost.get();
  List cost_optim_result;
  if (cost_gradient.isNull() && cost_hessian.isNull()) {
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
    cost_optim_result = List::create(
      Named("par") = std::log(
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

mat Fastcpd::get_theta_sum() {
  return theta_sum;
}

void Fastcpd::create_cost_function_wrapper(Nullable<Function> cost) {
  DEBUG_RCOUT(family);
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

void Fastcpd::create_cost_gradient_wrapper(Nullable<Function> cost_gradient) {
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

void Fastcpd::create_cost_hessian_wrapper(Nullable<Function> cost_hessian) {
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
    theta_hat.col(0) = segment_theta_hat.row(0).t();
    theta_sum.col(0) = segment_theta_hat.row(0).t();
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

void Fastcpd::create_segment_indices() {
  const unsigned int segment_length = floor(n / segment_count);
  const int segment_remainder = n % segment_count;
  for (
    int segment_index = 1; segment_index <= segment_count; segment_index++
  ) {
    if (segment_index <= segment_remainder) {
      segment_indices(span(
        (segment_index - 1) * (segment_length + 1),
        segment_index * (segment_length + 1)
      )).fill(segment_index);
    } else {
      segment_indices(span(
        (segment_index - 1) * segment_length + segment_remainder,
        segment_index * segment_length + segment_remainder - 1
      )).fill(segment_index);
    }
  }
}

// TODO(doccstat): Use `segment_theta` as warm start.

void Fastcpd::create_segment_statistics() {
  for (
    int segment_index = 0; segment_index < segment_count; ++segment_index
  ) {
    DEBUG_RCOUT(segment_index);
    ucolvec segment_indices_ = find(segment_indices == segment_index + 1);
    mat data_segment = data.rows(segment_indices_);
    rowvec segment_theta;
    if (!contain(FASTCPD_FAMILIES, family)) {
      segment_theta = as<rowvec>(get_optimized_cost(data_segment)["par"]);
    } else {
      segment_theta = as<rowvec>(
        cost_function_wrapper(
          data_segment, R_NilValue, 0, true, R_NilValue
        )["par"]
      );
    }
    DEBUG_RCOUT(segment_theta);

    // Initialize the estimated coefficients for each segment to be the
    // estimated coefficients in the segment.
    segment_theta_hat.row(segment_index) = segment_theta;
    if (family == "lasso" || family == "gaussian") {
      colvec segment_residual = data_segment.col(0) -
        data_segment.cols(1, data_segment.n_cols - 1) * segment_theta.t();
        double err_var = as_scalar(mean(square(segment_residual)));
        update_err_sd(segment_index, err_var);
        DEBUG_RCOUT(err_sd);
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

void Fastcpd::update_err_sd(
  const unsigned int segment_index, const double err_var
) {
  err_sd(segment_index) = sqrt(err_var);
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

void Fastcpd::update_theta_hat(colvec new_theta_hat) {
  theta_hat = join_rows(theta_hat, new_theta_hat);
}

void Fastcpd::update_theta_hat(
  const unsigned int col, colvec new_theta_hat
) {
  theta_hat.col(col) = new_theta_hat;
}

void Fastcpd::update_theta_sum(
  const unsigned int col, colvec new_theta_sum
) {
  theta_sum.col(col) += new_theta_sum;
}

void Fastcpd::update_theta_sum(colvec new_theta_sum) {
  theta_sum = join_rows(theta_sum, new_theta_sum);
}

void Fastcpd::update_theta_hat(ucolvec pruned_left) {
  theta_hat = theta_hat.cols(pruned_left);
}

void Fastcpd::update_theta_sum(ucolvec pruned_left) {
  theta_sum = theta_sum.cols(pruned_left);
}

List Fastcpd::run() {
  // Set up the initial values.
  double lambda = 0;
  mat start = zeros<mat>(p, n);

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

  if (vanilla_percentage < 1) {
    create_gradients();
  }

  checkUserInterrupt();
  if (r_progress) {
    rProgress.tick();
  }

#ifdef EXPERIMENT

  std::vector<double> cmat_vec((n + 1) * n, 0.0);
  colvec tau_stars = zeros<vec>(n + 1);

  std::for_each(
    std::execution::par_unseq,
    cmat_vec.begin(),
    cmat_vec.end(),
    [&](double& cval) {
      int index = &cval - &cmat_vec[0];
      int s = index % (n + 1);
      int t = index / (n + 1);
      if (s > t) {
        List cost_optim_result = negative_log_likelihood(
          data.rows(t, s - 1), R_NilValue, lambda, false, R_NilValue
        );
        cval = as<double>(cost_optim_result["value"]);
      } else {
        cval = arma::datum::inf;
      }
    }
  );

  arma::Mat<double> cmat(cmat_vec.data(), n + 1, n, /* copy_aux_mem */ false);

  for (int t = 1; t <= n; t++) {
    colvec new_fvec = fvec(t - 1) + cmat.col(t - 1) + beta;
    ucolvec f_t_condition = find(new_fvec < fvec);
    if (f_t_condition.n_elem > 0) {
      fvec.rows(f_t_condition) = new_fvec.rows(f_t_condition);
      tau_stars.rows(f_t_condition) = (t - 1) * ones<vec>(f_t_condition.n_elem);
    }
    cp_sets[t] = join_cols(cp_sets[tau_stars(t - 1)], colvec{tau_stars(t - 1)});
  }

#else

  for (int t = 1; t <= n; t++) {
    DEBUG_RCOUT(t);
    unsigned int r_t_count = r_t_set.n_elem;
    DEBUG_RCOUT(r_t_count);

    // Number of cost values is the same as the number of elements in R_t.
    colvec cval = zeros<vec>(r_t_count);

    // For tau in R_t \ {t-1}.
    for (unsigned int i = 1; i < r_t_count; i++) {
      DEBUG_RCOUT(i);
      int tau = r_t_set(i - 1);
      if (family == "lasso") {
        // Mean of `err_sd` only works if error sd is unchanged.
        lambda = mean(err_sd) * sqrt(2 * std::log(p) / (t - tau));
      }
      mat data_segment = data.rows(tau, t - 1);
      DEBUG_RCOUT(data_segment);
      if (t > vanilla_percentage * n) {
        // fastcpd
        update_cost_parameters(t, tau, i, k.get(), lambda, line_search);
        colvec theta = theta_sum.col(i - 1) / (t - tau);
        DEBUG_RCOUT(theta);
        if (!contain(FASTCPD_FAMILIES, family)) {
          Function cost_non_null = cost.get();
          SEXP cost_result = cost_non_null(data_segment, theta);
          cval(i - 1) = as<double>(cost_result);
        } else if (
          (family != "lasso" && t - tau >= p) ||
          (family == "lasso" && t - tau >= 3)
        ) {
          List cost_result = cost_function_wrapper(
            data_segment, wrap(theta), lambda, false, R_NilValue
          );
          cval(i - 1) = as<double>(cost_result["value"]);
        } else {
          // t - tau < p or for lasso t - tau < 3
        }
      } else {
        // vanilla PELT
        List cost_optim_result;
        if (!contain(FASTCPD_FAMILIES, family)) {
          cost_optim_result = get_optimized_cost(data_segment);
        } else {
          if (warm_start && t - tau >= 10 * p) {
            cost_optim_result =
              cost_function_wrapper(
                data_segment, R_NilValue, lambda, false,
                wrap(segment_theta_hat[segment_indices(t - 1) - 1])
                // Or use `wrap(start.col(tau))` for warm start.
            );
            start.col(tau) = as<colvec>(cost_optim_result["par"]);
          } else {
            cost_optim_result = cost_function_wrapper(
              data_segment, R_NilValue, lambda, false, R_NilValue
            );
          }
        }
        cval(i - 1) = as<double>(cost_optim_result["value"]);

        // If `vanilla_percentage` is not 1, then we need to keep track of
        // thetas for later `fastcpd` steps.
        if (vanilla_percentage < 1 && t <= vanilla_percentage * n) {
          update_theta_hat(i - 1, as<colvec>(cost_optim_result["par"]));
          update_theta_sum(i - 1, as<colvec>(cost_optim_result["par"]));
        }
      }
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
    ucolvec pruned_left = pruning ?
      find(cval + fvec.rows(r_t_set) + convexity_coef <= min_obj) : r_t_set;
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

#endif  // EXPERIMENT

  // Remove change points close to the boundaries.
  colvec raw_cp_set = cp_sets[n],
         cp_set = cp_sets[n];
  cp_set = cp_set(find(cp_set > trim * n));
  cp_set = cp_set(find(cp_set < (1 - trim) * n));
  colvec cp_set_ = zeros<vec>(cp_set.n_elem + 1);
  if (cp_set.n_elem) {
    cp_set_.rows(1, cp_set_.n_elem - 1) = std::move(cp_set);
  }
  cp_set = sort(unique(std::move(cp_set_)));

  // Remove change points close to each other.
  ucolvec cp_set_too_close = find(diff(cp_set) <= trim * n);
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
  cp_set = cp_set(find(cp_set > 0));

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
  cp_loc_(cp_loc_.n_elem - 1) = n;
  colvec cp_loc = unique(std::move(cp_loc_));
  colvec cost_values = zeros<vec>(cp_loc.n_elem - 1);
  mat thetas = zeros<mat>(p, cp_loc.n_elem - 1);
  mat residual;
  if (
    family == "mean" || family == "variance" ||
    family == "meanvariance" || family == "mv"
  ) {
    residual = zeros<mat>(data.n_rows, data.n_cols);
  } else if (family == "mgaussian") {
    residual = zeros<mat>(data.n_rows, p_response);
  } else {
    residual = zeros<mat>(data.n_rows, 1);
  }
  unsigned int residual_next_start = 0;

  for (unsigned int i = 0; i < cp_loc.n_elem - 1; i++) {
    colvec segment_data_index_ =
        linspace(cp_loc(i), cp_loc(i + 1) - 1, cp_loc(i + 1) - cp_loc(i));
    ucolvec segment_data_index =
        arma::conv_to<ucolvec>::from(std::move(segment_data_index_));

    mat data_segment = data.rows(segment_data_index);
    List cost_optim_result;
    if (!contain(FASTCPD_FAMILIES, family)) {
      cost_optim_result = get_optimized_cost(data_segment);
    } else {
      cost_optim_result = cost_function_wrapper(
        data_segment, R_NilValue, lambda, false, R_NilValue
      );
    }

    cost_values(i) = as<double>(cost_optim_result["value"]);

    // Parameters are not involved for PELT.
    if (vanilla_percentage < 1) {
      thetas.col(i) = as<colvec>(cost_optim_result["par"]);
    }

    // Residual is only calculated for built-in families.
    if (contain(FASTCPD_FAMILIES, family)) {
      mat cost_optim_residual = as<mat>(cost_optim_result["residuals"]);
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

void Fastcpd::update_cost_parameters(
  const unsigned int t,
  const unsigned int tau,
  const unsigned int i,
  Function k,
  const double lambda,
  const colvec line_search
) {
  List cost_update_result = update_cost_parameters_steps(
    data.rows(0, t - 1), tau, i, k, momentum, lambda, line_search
  );
  DEBUG_RCOUT(line_search);
  update_theta_hat(i - 1, as<colvec>(cost_update_result[0]));
  create_theta_sum(i - 1, as<colvec>(cost_update_result[1]));
  update_hessian(i - 1, as<mat>(cost_update_result[2]));
  update_momentum(as<colvec>(cost_update_result[3]));
}

List Fastcpd::update_cost_parameters_steps(
    const mat data,
    const int tau,
    const int i,
    Function k,
    colvec momentum,
    const double lambda,
    colvec line_search
) {
  update_cost_parameters_step(data, i, 0, data.n_rows - 1, lambda, line_search);

  for (int kk = 1; kk <= as<int>(k(data.n_rows - tau)); kk++) {
    for (unsigned j = tau + 1; j <= data.n_rows; j++) {
      update_cost_parameters_step(data, i, tau, j - 1, lambda, line_search);
    }
  }

  theta_sum.col(i - 1) += theta_hat.col(i - 1);
  return List::create(
    theta_hat.col(i - 1), theta_sum.col(i - 1), hessian.slice(i - 1), momentum
  );
}

void Fastcpd::update_cost_parameters_step(
  const mat data,
  const int i,
  const int data_start,
  const int data_end,
  const double lambda,
  const colvec line_search
) {
  DEBUG_RCOUT(data_start);
  mat hessian_i = hessian.slice(i - 1);
  colvec gradient;

  if (!contain(FASTCPD_FAMILIES, family)) {
    mat cost_hessian_result = cost_hessian_wrapper(
      data.rows(data_start, data_end), theta_hat.col(i - 1)
    );
    DEBUG_RCOUT(cost_hessian_result);
    hessian_i += cost_hessian_result;
    colvec cost_gradient_result = cost_gradient_wrapper(
      data.rows(data_start, data_end), theta_hat.col(i - 1)
    );
    gradient = cost_gradient_result;
    DEBUG_RCOUT(gradient);
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
  DEBUG_RCOUT(momentum);

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
        theta_hat.col(i - 1) + line_search[line_search_index] * momentum;
      DEBUG_RCOUT(theta_candidate);
      colvec theta_upper_bound = arma::min(std::move(theta_candidate), upper);
      colvec theta_projected = arma::max(std::move(theta_upper_bound), lower);
      line_search_costs[line_search_index] = cost_function_wrapper(
        data, wrap(theta_projected), lambda, false, R_NilValue
      )["value"];
    }
  }
  best_learning_rate = line_search[line_search_costs.index_min()];
  DEBUG_RCOUT(best_learning_rate);

  // Update theta_hat with momentum
  theta_hat.col(i - 1) += best_learning_rate * momentum;

  theta_hat.col(i - 1) = arma::min(theta_hat.col(i - 1), upper);
  theta_hat.col(i - 1) = arma::max(theta_hat.col(i - 1), lower);

  if (family == "lasso" || family == "gaussian") {
    // Update theta_hat with L1 penalty
    double hessian_norm = norm(hessian_i, "fro");
    vec normd = abs(theta_hat.col(i - 1)) - lambda / hessian_norm;
    theta_hat.col(i - 1) = sign(theta_hat.col(i - 1)) % arma::max(
      normd, zeros<colvec>(normd.n_elem)
    );
  }

  hessian.slice(i - 1) = std::move(hessian_i);
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
    cum_coef_add = coef_add;
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

void Fastcpd::update_momentum(colvec new_momentum) {
  momentum = new_momentum;
}

}  // namespace fastcpd::classes
