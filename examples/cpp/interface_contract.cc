#include <fastcpd/fastcpd.h>

#include <armadillo>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;

fs::path find_fixture_directory(int argc, char** argv) {
  std::vector<fs::path> candidates;
  if (argc > 1) candidates.emplace_back(argv[1]);
  candidates.emplace_back("tests/fixtures");
  candidates.emplace_back("../tests/fixtures");
  candidates.emplace_back("../../tests/fixtures");
  for (fs::path const& candidate : candidates) {
    std::error_code error;
    if (fs::is_regular_file(candidate / "manifest.tsv", error)) {
      return candidate;
    }
  }
  throw std::runtime_error(
      "unable to locate tests/fixtures; pass its path as the first argument");
}

arma::mat read_csv(fs::path const& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("unable to open " + path.string());
  std::string line;
  if (!std::getline(input, line)) {
    throw std::runtime_error("empty CSV fixture: " + path.string());
  }
  std::size_t columns =
      static_cast<std::size_t>(std::count(line.begin(), line.end(), ',')) + 1;
  std::vector<double> values;
  std::size_t rows = 0;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::stringstream stream(line);
    std::string field;
    std::size_t row_columns = 0;
    while (std::getline(stream, field, ',')) {
      values.push_back(std::stod(field));
      ++row_columns;
    }
    if (row_columns != columns) {
      throw std::runtime_error("inconsistent CSV fixture: " + path.string());
    }
    ++rows;
  }
  arma::mat result(rows, columns);
  for (std::size_t row = 0; row < rows; ++row) {
    for (std::size_t column = 0; column < columns; ++column) {
      result(row, column) = values[row * columns + column];
    }
  }
  return result;
}

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

int main(int argc, char** argv) {
  try {
    fs::path const fixtures = find_fixture_directory(argc, argv);
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

    // Omitted time-series orders use the same defaults in generic and named
    // entry points. Explicit all-zero orders remain invalid.
    arma::colvec const default_series =
        arma::sin(arma::linspace<arma::colvec>(0.0, 23.0, 24));
    std::vector<std::pair<std::string, arma::colvec>> const default_orders = {
        {"ar", {1.0}}, {"arma", {1.0, 1.0}},
        {"arima", {1.0, 1.0, 0.0}}, {"garch", {1.0, 1.0}}};
    for (auto const& entry : default_orders) {
      fastcpd::Options options = mean_options();
      options.beta = 1e6;
      options.family = entry.first;
      fastcpd::Result generic = fastcpd::detect(default_series, options);
      options.family = "mean";  // Named wrappers choose their own family.
      fastcpd::Result named = entry.first == "ar"
          ? fastcpd::detect_ar(default_series, options)
          : entry.first == "arma"
              ? fastcpd::detect_arma(default_series, options)
              : entry.first == "arima"
                  ? fastcpd::detect_arima(default_series, options)
                  : fastcpd::detect_garch(default_series, options);
      if (!arma::approx_equal(generic.order, entry.second, "absdiff", 0.0) ||
          !arma::approx_equal(named.order, entry.second, "absdiff", 0.0) ||
          !arma::approx_equal(generic.change_points, named.change_points,
                              "absdiff", 0.0)) {
        throw std::runtime_error("default order contract failed: " + entry.first);
      }
      options.family = entry.first;
      options.order.zeros(entry.second.n_elem);
      expect_invalid("explicit zero " + entry.first + " order", [&] {
        fastcpd::detect(default_series, options);
      });
    }

    arma::mat default_design(default_series.n_rows, 2, arma::fill::ones);
    default_design.col(0) = default_series;
    fastcpd::Options default_options;
    default_options.beta = 1e6;
    for (std::string const family : {"lm", "lasso", "binomial", "poisson",
                                     "quantile", "rank"}) {
      arma::mat family_data = default_design;
      if (family == "binomial" || family == "poisson") {
        for (arma::uword row = 0; row < family_data.n_rows; ++row) {
          family_data(row, 0) = static_cast<double>(row % 2);
        }
      }
      if (family == "rank") family_data = default_series;
      default_options.family = family;
      fastcpd::Result const result = fastcpd::detect(family_data, default_options);
      arma::colvec const expected_order = family == "quantile"
          ? arma::colvec{0.5} : arma::colvec{0.0, 0.0, 0.0};
      if (!arma::approx_equal(result.order, expected_order, "absdiff", 0.0)) {
        throw std::runtime_error("default metadata contract failed: " + family);
      }
    }
    default_options.family = "mean";
    fastcpd::Result const default_quantile =
        fastcpd::detect_quantile(default_design, default_options);
    if (default_quantile.order.n_elem != 1 || default_quantile.order(0) != 0.5) {
      throw std::runtime_error("named quantile did not default to the median");
    }
    arma::mat default_var_data(default_series.n_rows, 2);
    default_var_data.col(0) = default_series;
    default_var_data.col(1) = arma::cos(
        arma::linspace<arma::colvec>(0.0, 23.0, 24) * 0.7);
    fastcpd::Result const default_var =
        fastcpd::detect_var(default_var_data, default_options);
    default_options.family = "var";
    fastcpd::Result const generic_var = fastcpd::detect(default_var_data, default_options);
    if (default_var.order.n_elem != 1 || default_var.order(0) != 1.0 ||
        !arma::approx_equal(default_var.order, generic_var.order, "absdiff", 0.0)) {
      throw std::runtime_error("default VAR order contract failed");
    }
    default_options.family = "kcp";
    default_options.seed = 31;
    fastcpd::Result const default_kcp = fastcpd::detect(default_series, default_options);
    if (!arma::approx_equal(default_kcp.order, arma::colvec{100.0, 0.0},
                            "absdiff", 0.0)) {
      throw std::runtime_error("default KCP order contract failed");
    }

    fastcpd::ConfidenceOptions empty_profile_options;
    empty_profile_options.method = "profile";
    empty_profile_options.min_segment_length = 5;
    empty_profile_options.detector_options = mean_options();
    std::vector<fastcpd::ConfidenceInterval> const empty_profile =
        fastcpd::confint(named_mean, data, empty_profile_options);
    if (empty_profile.size() != 1 ||
        !std::isnan(empty_profile.front().lower) ||
        !std::isnan(empty_profile.front().upper)) {
      throw std::runtime_error("short profile candidate range was not empty");
    }

    expect_invalid("confidence result boundary", [&] {
      fastcpd::Result invalid = named_mean;
      invalid.change_points(0) = static_cast<double>(data.n_rows);
      fastcpd::confint(invalid, data, empty_profile_options);
    });

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

    arma::mat kernel_data(24, 2);
    for (arma::uword row = 0; row < kernel_data.n_rows; ++row) {
      kernel_data(row, 0) = static_cast<double>(row + 1) / 10.0;
      kernel_data(row, 1) = row % 2 == 0 ? -1.0 : 1.0;
    }
    fastcpd::Options kernel_options = mean_options();
    kernel_options.order = arma::colvec{8.0, 1.25};
    kernel_options.seed = 7;
    kernel_options.beta = 2.0;
    kernel_options.cp_only = false;
    kernel_options.variance_estimate.reset();
    fastcpd::Result kernel =
        fastcpd::detect_kernel(kernel_data, kernel_options);
    if (kernel.family != "kcp" || kernel.change_points.n_elem != 1 ||
        kernel.change_points(0) != 13.0 || kernel.cost_values.n_elem != 2 ||
        std::abs(kernel.cost_values(0) - 2.6735583651891526) > 1e-10 ||
        std::abs(kernel.cost_values(1) - 1.9910581518133954) > 1e-10) {
      throw std::runtime_error("seeded KCP contract failed");
    }

    fastcpd::Options constant_kernel_options = mean_options();
    constant_kernel_options.order = arma::colvec{4.0, 0.0};
    constant_kernel_options.seed = 7;
    constant_kernel_options.variance_estimate.reset();
    fastcpd::Result constant_kernel = fastcpd::detect_kernel(
        arma::zeros<arma::mat>(12, 1), constant_kernel_options);
    if (!constant_kernel.cost_values.is_finite()) {
      throw std::runtime_error("constant-input KCP bandwidth fallback failed");
    }

    arma::mat const mean_fixture = read_csv(fixtures / "mean_step.csv");
    arma::mat const mean_variance =
        fastcpd::estimate_variance_mean(mean_fixture);
    if (mean_variance.n_rows != 1 ||
        std::abs(mean_variance(0, 0) - 0.3205128205128205) > 1e-12) {
      throw std::runtime_error("mean variance fixture failed");
    }

    arma::mat const rank_fixture = read_csv(fixtures / "rank_step.csv");
    double const median_variance =
        fastcpd::estimate_variance_median(rank_fixture);
    if (std::abs(median_variance - 13.32515158156184) > 1e-12) {
      throw std::runtime_error("median variance fixture failed");
    }

    arma::mat const lm_fixture = read_csv(fixtures / "lm_step.csv");
    arma::mat const lm_variance =
        fastcpd::estimate_variance_linear_regression(lm_fixture);
    if (lm_variance.n_rows != 1 ||
        std::abs(lm_variance(0, 0) - 97.95856069973233) > 1e-9) {
      throw std::runtime_error("linear-regression variance fixture failed");
    }

    arma::colvec const arma_fixture =
        read_csv(fixtures / "arma_variance.csv").col(0);
    fastcpd::VarianceArmaResult const arma_variance =
        fastcpd::estimate_variance_arma(arma_fixture, 2, 2);
    if (arma_variance.table.size() != 4 ||
        std::abs(arma_variance.sigma2_aic - 0.8938872868847426) > 2e-6 ||
        std::abs(arma_variance.sigma2_bic - 0.8938872868847426) > 2e-6) {
      throw std::runtime_error("ARMA variance fixture failed");
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
