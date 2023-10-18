#include "constants.h"
#include "cost_update.h"
#include "functions.h"

using ::fastcpd::functions::cost_update_gradient;
using ::fastcpd::functions::cost_update_hessian;

namespace fastcpd::cost_update {

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

  if (contain(CUSTOM_FAMILIES, family)) {
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
      if (contain(CUSTOM_FAMILIES, family)) {
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
  List cost_optim_result;
  if (family == "vanilla") {
    cost_optim_result = List::create(
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
      Named("fn") = cost,
      Named("method") = "L-BFGS-B",
      Named("data") = data_segment
    );
    cost_optim_result = List::create(
      Named("par") = optim_result["par"],
      Named("value") = optim_result["value"],
      Named("residuals") = R_NilValue
    );
  } else {
    // This branch should not be reached.
  }
  return cost_optim_result;
}

}  // namespace fastcpd::cost_update
