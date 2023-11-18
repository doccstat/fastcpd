#include "constants.h"
#include "cost_update.h"
#include "functions.h"

using ::fastcpd::functions::cost_update_gradient;
using ::fastcpd::functions::cost_update_hessian;

namespace fastcpd::cost_update {

List cost_optim(
    const string family,
    const double vanilla_percentage,
    const int p,
    const mat data_segment,
    Function cost,
    const double lambda,
    const bool cv,
    const colvec lower,
    const colvec upper
) {
  List cost_optim_result;
  if (vanilla_percentage == 1) {
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
    stop("This branch should not be reached at cost_update.cc: 198.");
    // # nocov end
  }
  return cost_optim_result;
}

}  // namespace fastcpd::cost_update
