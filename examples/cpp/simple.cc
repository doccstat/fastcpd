#include <fastcpd/fastcpd.h>

#include <armadillo>

#include <cmath>
#include <iostream>

int main() {
  arma::mat data(100, 1);
  data.rows(0, 49).fill(0.0);
  data.rows(50, 99).fill(5.0);

  fastcpd::Options options;
  options.family = "mean";
  options.beta = 5.0;
  options.cost_adjustment = "BIC";
  options.cp_only = true;
  options.variance_estimate = arma::eye(1, 1);

  fastcpd::Result result = fastcpd::detect(data, options);

  result.change_points.t().print(std::cout, "change points:");
  if (result.change_points.n_elem != 1 ||
      std::abs(result.change_points(0) - 50.0) > 1.0) {
    return 1;
  }
  return 0;
}
