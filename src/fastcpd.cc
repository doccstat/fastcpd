#include "constants.h"
#include "fastcpd.h"
#include "functions.h"
#include "parameters.h"
#include "wrappers.h"

using ::fastcpd::functions::cost_update_gradient;
using ::fastcpd::functions::cost_update_hessian;
using ::fastcpd::functions::negative_log_likelihood;

List cost_update(
    const mat data,
    mat theta_hat,
    mat theta_sum,
    cube hessian,
    const int tau,
    const int i,
    Function k,
    const string family,
    colvec momentum,
    const double momentum_coef,
    const double epsilon,
    const double min_prob,
    const double winsorise_minval,
    const double winsorise_maxval,
    const double lambda,
    function<colvec(mat data, colvec theta, string family)>
      cost_gradient_wrapper,
    function<mat(mat data, colvec theta, string family, double min_prob)>
      cost_hessian_wrapper
) {
  // Get the hessian
  mat hessian_i = hessian.slice(i - 1);
  colvec gradient;

  if (CUSTOM_FAMILIES.find(family) != CUSTOM_FAMILIES.end()) {
    mat cost_hessian_result = cost_hessian_wrapper(
      data, theta_hat.col(i - 1),
      family,  // UNUSED
      min_prob  // UNUSED
    );
    colvec cost_gradient_result = cost_gradient_wrapper(
      data, theta_hat.col(i - 1),
      family  // UNUSED
    );
    hessian_i += cost_hessian_result;
    gradient = cost_gradient_result;
  } else {
    hessian_i +=
      cost_update_hessian(data, theta_hat.col(i - 1), family, min_prob);
    gradient = cost_update_gradient(data, theta_hat.col(i - 1), family);
  }

  // Add epsilon to the diagonal for PSD hessian
  mat hessian_psd = hessian_i + epsilon * eye<mat>(
    theta_hat.n_rows, theta_hat.n_rows
  );

  // Calculate momentum step
  vec momentum_step = solve(hessian_psd, gradient);
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
    theta_hat.col(i - 1) =
      vec(winsorize_result.begin(), winsorize_result.size(), false);
  } else if (family == "lasso" || family == "gaussian") {
    // Update theta_hat with L1 penalty
    double hessian_norm = norm(hessian_i, "fro");
    vec normd = arma::abs(theta_hat.col(i - 1)) - lambda / hessian_norm;
    theta_hat.col(i - 1) =
        sign(theta_hat.col(i - 1)) % max(normd, zeros<colvec>(normd.n_elem));
  }

  for (int kk = 1; kk <= as<int>(k(data.n_rows - tau)); kk++) {
    for (unsigned j = tau + 1; j <= data.n_rows; j++) {
      if (CUSTOM_FAMILIES.find(family) != CUSTOM_FAMILIES.end()) {
        mat cost_hessian_result = cost_hessian_wrapper(
          data.rows(tau, j - 1), theta_hat.col(i - 1),
          family,  // UNUSED
          min_prob  // UNUSED
        );
        hessian_i += cost_hessian_result;
        colvec cost_gradient_result = cost_gradient_wrapper(
          data.rows(tau, j - 1), theta_hat.col(i - 1),
          family  // UNUSED
        );
        gradient = cost_gradient_result;
      } else {
        hessian_i += cost_update_hessian(
          data.rows(tau, j - 1), theta_hat.col(i - 1), family, min_prob
        );
        gradient = cost_update_gradient(
          data.rows(tau, j - 1), theta_hat.col(i - 1), family
        );
      }

      hessian_psd =
        hessian_i + epsilon * eye<mat>(theta_hat.n_rows, theta_hat.n_rows);
      momentum_step = solve(hessian_psd, gradient);
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
            vec(winsorize_result.begin(), winsorize_result.size(), false);
      } else if (family == "lasso" || family == "gaussian") {
        double hessian_norm = norm(hessian_i, "fro");
        vec normd = arma::abs(theta_hat.col(i - 1)) - lambda / hessian_norm;
        theta_hat.col(i - 1) =
          sign(theta_hat.col(i - 1)) % max(normd, zeros<colvec>(normd.n_elem));
      }
    }
  }

  theta_sum.col(i - 1) += theta_hat.col(i - 1);
  hessian.slice(i - 1) = std::move(hessian_i);
  return List::create(
    theta_hat.col(i - 1), theta_sum.col(i - 1), hessian.slice(i - 1), momentum
  );
}

List cost_optim(
    const string family,
    const int p,
    const mat data_segment,
    Function cost,
    const double lambda,
    const bool cv
) {
  if (family == "vanilla") {
    return List::create(
      Named("par") = R_NilValue,
      Named("value") = cost(data_segment),
      Named("residuals") = R_NilValue
    );
  } else if (p == 1) {
    Environment stats = Environment::namespace_env("stats");
    Function optim = stats["optim"];
    List optim_result = optim(
      Named("par") = 0,
      Named("fn") = InternalFunction(
        +[](double theta, mat data, Function cost) {
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
  } else if (p > 1) {
    Environment stats = Environment::namespace_env("stats");
    Function optim = stats["optim"];
    List optim_result = optim(
      Named("par") = zeros<vec>(p),
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
    // NOTEST
    stop("Internal error, this branch should not be reached.");
  }
}

List fastcpd_impl(
    mat data,
    double beta,
    const int segment_count,
    const double trim,
    const double momentum_coef,
    Function k,
    const string family,
    const double epsilon,
    const double min_prob,
    const double winsorise_minval,
    const double winsorise_maxval,
    const int p,
    Nullable<Function> cost,
    Nullable<Function> cost_gradient,
    Nullable<Function> cost_hessian,
    const bool cp_only,
    const double vanilla_percentage
) {
  // Set up the initial values.
  const int n = data.n_rows;
  double lambda = 0;

  // After t = 1, the r_t_set R_t contains 0 and 1.
  ucolvec r_t_set = {0, 1};
  // C(0) = NULL, C(1) = {0}.
  std::vector<colvec> cp_sets = {zeros<vec>(0)};
  linspace(1, n, n).for_each([&](int i) {
    cp_sets.push_back(zeros<vec>(1));
  });
  // Objective function: F(0) = -beta.
  colvec f_t = zeros<vec>(n + 1);
  f_t(0) = -beta;

  fastcpd::parameters::FastcpdParameters fastcpd_parameters_class(
    data, beta, p, family, segment_count,
    winsorise_minval, winsorise_maxval, epsilon
  );

  fastcpd_parameters_class.wrap_cost(cost);
  fastcpd_parameters_class.wrap_cost_gradient(cost_gradient);
  fastcpd_parameters_class.wrap_cost_hessian(cost_hessian);

  if (vanilla_percentage != 1) {
    fastcpd_parameters_class.create_segment_indices();
    fastcpd_parameters_class.create_segment_statistics();
    fastcpd_parameters_class.update_beta();
    fastcpd_parameters_class.create_gradients();
  }

  for (int t = 2; t <= n; t++) {
    unsigned int r_t_count = r_t_set.n_elem;

    // Number of cost values is the same as the number of elements in R_t.
    colvec cval = zeros<vec>(r_t_count);

    // For tau in R_t \ {t-1}.
    for (unsigned int i = 1; i < r_t_count; i++) {
      int tau = r_t_set(i - 1);
      if (family == "lasso") {
        // Mean of `err_sd` only works if error sd is unchanged.
        lambda = mean(
          fastcpd_parameters_class.get_err_sd()
        ) * sqrt(2 * log(p) / (t - tau));
      }
      mat data_segment = data.rows(tau, t - 1);
      if (vanilla_percentage == 1 || t <= vanilla_percentage * n) {
        List cost_optim_result;
        if (CUSTOM_FAMILIES.find(family) != CUSTOM_FAMILIES.end()) {
          cost_optim_result = cost_optim(
            family, p, data_segment,
            fastcpd_parameters_class.cost.get(), 0, true
          );
        } else {
          cost_optim_result =
            negative_log_likelihood(data_segment, R_NilValue, family, 0, true);
        }
        cval(i - 1) = as<double>(cost_optim_result["value"]);
        if (vanilla_percentage < 1 && t <= vanilla_percentage * n) {
          fastcpd_parameters_class.update_theta_hat(
            i - 1, as<colvec>(cost_optim_result["par"])
          );
          fastcpd_parameters_class.update_theta_sum(
            i - 1, as<colvec>(cost_optim_result["par"])
          );
        }
      } else {
        List cost_update_result = cost_update(
          data.rows(0, t - 1),
          fastcpd_parameters_class.get_theta_hat(),
          fastcpd_parameters_class.get_theta_sum(),
          fastcpd_parameters_class.get_hessian(),
          tau,
          i,
          k,
          family,
          fastcpd_parameters_class.get_momentum(),
          momentum_coef,
          epsilon,
          min_prob,
          winsorise_minval,
          winsorise_maxval,
          lambda,
          fastcpd_parameters_class.cost_gradient_wrapper,
          fastcpd_parameters_class.cost_hessian_wrapper
        );
        fastcpd_parameters_class.update_theta_hat(
          i - 1, as<colvec>(cost_update_result[0])
        );
        fastcpd_parameters_class.create_theta_sum(
          i - 1, as<colvec>(cost_update_result[1])
        );
        fastcpd_parameters_class.update_hessian(
          i - 1, as<mat>(cost_update_result[2])
        );
        fastcpd_parameters_class.update_momentum(
          as<colvec>(cost_update_result[3])
        );
        colvec theta =
            fastcpd_parameters_class.get_theta_sum().col(i - 1) / (t - tau);
        if (family == "poisson" && t - tau >= p) {
          Environment desc_tools = Environment::namespace_env("DescTools");
          Function winsorize = desc_tools["Winsorize"];
          NumericVector winsorize_result = winsorize(
            Rcpp::_["x"] = theta,
            Rcpp::_["minval"] = winsorise_minval,
            Rcpp::_["maxval"] = winsorise_maxval
          );
          theta = as<colvec>(winsorize_result);
        }
        if (
          (
            FASTCPD_FAMILIES.find(family) != FASTCPD_FAMILIES.end() &&
            family != "lasso" && t - tau >= p
          ) ||
          (family == "lasso" && t - tau >= 3)
        ) {
          List cost_result = negative_log_likelihood(
            data_segment, Rcpp::wrap(theta), family, lambda
          );
          cval(i - 1) = as<double>(cost_result["value"]);
        } else if (CUSTOM_FAMILIES.find(family) != CUSTOM_FAMILIES.end()) {
          // if (warm_start && t - tau >= 50) {
          //   cost_result <- cost(data_segment, start = start[, tau + 1])
          //   start[, tau + 1] <- cost_result$par
          //   cval[i] <- cost_result$value
          // } else {
          Function cost_non_null = fastcpd_parameters_class.cost.get();
          SEXP cost_result = cost_non_null(data_segment, theta);
          cval(i - 1) = as<double>(cost_result);
          // }
        }
      }
    }

    if (vanilla_percentage != 1) {
      fastcpd_parameters_class.update_fastcpd_parameters(t);
    }

    // Step 3
    cval(r_t_count - 1) = 0;

    // `beta` adjustment seems to work but there might be better choices.
    colvec obj = cval + f_t.rows(r_t_set) + beta;
    double min_obj = min(obj);
    double tau_star = r_t_set(index_min(obj));

    // Step 4
    cp_sets[t] = join_cols(cp_sets[tau_star], colvec{tau_star});

    // Step 5
    ucolvec pruned_left = arma::find(cval + f_t.rows(r_t_set) <= min_obj);
    ucolvec pruned_r_t_set = zeros<ucolvec>(pruned_left.n_elem + 1);
    pruned_r_t_set.rows(0, pruned_left.n_elem - 1) = r_t_set(pruned_left);
    pruned_r_t_set(pruned_left.n_elem) = t;
    r_t_set = std::move(pruned_r_t_set);

    if (vanilla_percentage != 1) {
      fastcpd_parameters_class.update_theta_hat(pruned_left);
      fastcpd_parameters_class.update_theta_sum(pruned_left);
      fastcpd_parameters_class.update_hessian(pruned_left);
    }

    // Objective function F(t).
    f_t(t) = min_obj;
  }

  // Remove change points close to the boundaries.
  colvec cp_set = cp_sets[n];
  cp_set = cp_set(arma::find(cp_set >= trim * n));
  cp_set = cp_set(arma::find(cp_set <= (1 - trim) * n));
  colvec cp_set_ = zeros<vec>(cp_set.n_elem + 1);
  if (cp_set.n_elem) {
    cp_set_.rows(1, cp_set_.n_elem - 1) = std::move(cp_set);
  }
  cp_set = arma::sort(arma::unique(std::move(cp_set_)));

  // Remove change points close to each other.
  ucolvec cp_set_too_close = arma::find(arma::diff(cp_set) <= trim * n);
  if (cp_set_too_close.n_elem > 0) {
    int rest_element_count = cp_set.n_elem - cp_set_too_close.n_elem;
    colvec cp_set_rest_left = zeros<vec>(rest_element_count),
          cp_set_rest_right = zeros<vec>(rest_element_count);
    for (unsigned int i = 0, i_left = 0, i_right = 0; i < cp_set.n_elem; i++) {
      if (
        ucolvec left_find = arma::find(cp_set_too_close == i);
        left_find.n_elem == 0
      ) {
        cp_set_rest_left(i_left) = cp_set(i);
        i_left++;
      }
      if (
        ucolvec right_find = arma::find(cp_set_too_close == i - 1);
        right_find.n_elem == 0
      ) {
        cp_set_rest_right(i_right) = cp_set(i);
        i_right++;
      }
    }
    cp_set = arma::floor((cp_set_rest_left + cp_set_rest_right) / 2);
  }
  cp_set = cp_set(arma::find(cp_set > 0));

  if (cp_only) {
    return List::create(
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
  colvec cp_loc = arma::unique(std::move(cp_loc_));
  colvec cost_values = zeros<vec>(cp_loc.n_elem - 1);
  mat thetas = zeros<mat>(p, cp_loc.n_elem - 1);
  colvec residual = zeros<colvec>(n);
  unsigned int residual_next_start = 0;
  for (unsigned int i = 0; i < cp_loc.n_elem - 1; i++) {
    colvec segment_data_index_ =
        linspace(cp_loc(i), cp_loc(i + 1) - 1, cp_loc(i + 1) - cp_loc(i));
    ucolvec segment_data_index =
        arma::conv_to<ucolvec>::from(std::move(segment_data_index_));
    mat data_segment = data.rows(segment_data_index);
    List cost_optim_result;
    if (CUSTOM_FAMILIES.find(family) != CUSTOM_FAMILIES.end()) {
      cost_optim_result = cost_optim(
        family, p, data_segment,
        fastcpd_parameters_class.cost.get(), lambda, false
      );
    } else {
      cost_optim_result = negative_log_likelihood(
        data_segment, R_NilValue, family, lambda, false
      );
    }
    colvec cost_optim_par = as<colvec>(cost_optim_result["par"]);
    double cost_optim_value = as<double>(cost_optim_result["value"]);
    colvec cost_optim_residual = as<colvec>(cost_optim_result["residuals"]);
    thetas.col(i) = cost_optim_par;
    cost_values(i) = cost_optim_value;
    residual.rows(
      residual_next_start, residual_next_start + cost_optim_residual.n_elem - 1
    ) = cost_optim_residual;
    residual_next_start += cost_optim_residual.n_elem;
  }
  return List::create(
    Named("cp_set") = cp_set,
    Named("cost_values") = cost_values,
    Named("residual") = residual,
    Named("thetas") = thetas
  );
}
