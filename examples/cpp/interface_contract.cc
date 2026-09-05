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
    if (valid.family != "mean" || !valid.cp_only) {
      throw std::runtime_error("generic result metadata is incomplete");
    }

    fastcpd::Result named_mean = fastcpd::detect_mean(data, mean_options());
    if (named_mean.family != "mean" ||
        !arma::approx_equal(named_mean.change_points, valid.change_points,
                            "absdiff", 0.0)) {
      throw std::runtime_error("detect_mean does not match generic mean");
    }

    fastcpd::Options rank_options = mean_options();
    fastcpd::Result rank = fastcpd::detect_rank(data, rank_options);
    if (rank.family != "rank" || rank.change_points.n_elem != 1 ||
        rank.change_points(0) != 4.0) {
      throw std::runtime_error("rank wrapper contract failed");
    }

    arma::colvec series = arma::linspace<arma::colvec>(1.0, 16.0, 16);
    fastcpd::Options ar_options = mean_options();
    ar_options.order = arma::colvec{1.0};
    ar_options.beta = 1e6;
    ar_options.cp_only = false;
    fastcpd::Result ar = fastcpd::detect_ar(series, ar_options);
    if (ar.family != "ar" || ar.residuals.n_rows != series.n_rows ||
        !std::isnan(ar.residuals(0, 0))) {
      throw std::runtime_error("AR lag-coordinate restoration failed");
    }

    fastcpd::Options arma_options = ar_options;
    arma_options.order = arma::colvec{1.0, 0.0};
    fastcpd::Result pure_ar = fastcpd::detect_arma(series, arma_options);
    if (pure_ar.family != "arma" ||
        !arma::approx_equal(pure_ar.change_points, ar.change_points,
                            "absdiff", 0.0) ||
        pure_ar.residuals.n_rows != series.n_rows) {
      throw std::runtime_error("pure-AR ARMA routing failed");
    }

    arma::mat var_data(series.n_rows, 2);
    var_data.col(0) = series;
    var_data.col(1) = 2.0 * series + 1.0;
    fastcpd::Options var_options = mean_options();
    var_options.order = arma::colvec{1.0};
    var_options.beta = 1e6;
    var_options.cp_only = false;
    var_options.variance_estimate = arma::eye<arma::mat>(2, 2);
    fastcpd::Result var = fastcpd::detect_var(var_data, var_options);
    if (var.family != "var" || var.residuals.n_rows != var_data.n_rows ||
        var.residuals.n_cols != var_data.n_cols ||
        !std::isnan(var.residuals(0, 0))) {
      throw std::runtime_error("VAR raw-input contract failed");
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
