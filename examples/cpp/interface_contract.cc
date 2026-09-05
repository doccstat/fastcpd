#include <fastcpd/fastcpd.h>

#include <armadillo>

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

void expect_invalid(std::string const& label,
                    std::function<void()> const& operation) {
  try {
    operation();
  } catch (std::invalid_argument const&) {
    return;
  }
  throw std::runtime_error(label + " did not reject invalid input");
}

fastcpd::Options mean_options() {
  fastcpd::Options options;
  options.family = "mean";
  options.beta = 2.0;
  options.cost_adjustment = "BIC";
  options.cp_only = true;
  options.variance_estimate = arma::eye<arma::mat>(1, 1);
  return options;
}

}  // namespace

int main() {
  try {
    arma::mat data(8, 1, arma::fill::zeros);
    data.rows(4, 7).fill(3.0);

    fastcpd::Result valid = fastcpd::detect(data, mean_options());
    if (valid.change_points.n_elem != 1 || valid.change_points(0) != 4.0) {
      throw std::runtime_error("valid mean contract returned wrong change point");
    }

    expect_invalid("non-finite data", [&] {
      arma::mat invalid = data;
      invalid(0, 0) = std::numeric_limits<double>::quiet_NaN();
      fastcpd::detect(invalid, mean_options());
    });
    expect_invalid("trim", [&] {
      fastcpd::Options options = mean_options();
      options.trim = 1.1;
      fastcpd::detect(data, options);
    });
    expect_invalid("vanilla_percentage", [&] {
      fastcpd::Options options = mean_options();
      options.vanilla_percentage = -0.1;
      fastcpd::detect(data, options);
    });
    expect_invalid("segment_count", [&] {
      fastcpd::Options options = mean_options();
      options.segment_count = 0;
      fastcpd::detect(data, options);
    });
    expect_invalid("epsilon", [&] {
      fastcpd::Options options = mean_options();
      options.epsilon = 0.0;
      fastcpd::detect(data, options);
    });
    expect_invalid("numeric beta", [&] {
      fastcpd::Options options = mean_options();
      options.beta = std::numeric_limits<double>::infinity();
      fastcpd::detect(data, options);
    });
    expect_invalid("line_search", [&] {
      fastcpd::Options options = mean_options();
      options.line_search = arma::colvec{1.0, 0.0};
      fastcpd::detect(data, options);
    });
    expect_invalid("NaN lower", [&] {
      fastcpd::Options options = mean_options();
      options.lower = arma::colvec{
          std::numeric_limits<double>::quiet_NaN()};
      fastcpd::detect(data, options);
    });
    expect_invalid("inverted bounds", [&] {
      fastcpd::Options options = mean_options();
      options.lower = arma::colvec{2.0};
      options.upper = arma::colvec{1.0};
      fastcpd::detect(data, options);
    });
    expect_invalid("ARMA order", [&] {
      fastcpd::Options options = mean_options();
      options.family = "arma";
      options.order = arma::colvec{1.5, 1.0};
      fastcpd::detect(data, options);
    });
    expect_invalid("multivariate GARCH", [&] {
      fastcpd::Options options = mean_options();
      options.family = "garch";
      options.order = arma::colvec{1.0, 1.0};
      arma::mat multivariate = arma::join_rows(data, data);
      fastcpd::detect(multivariate, options);
    });
    expect_invalid("quantile level", [&] {
      fastcpd::Options options = mean_options();
      options.family = "quantile";
      options.order = arma::colvec{1.0};
      arma::mat regression =
          arma::join_rows(data, arma::ones<arma::mat>(8, 1));
      fastcpd::detect(regression, options);
    });

    std::cout << "standalone interface contract passed\n";
    return 0;
  } catch (std::exception const& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
