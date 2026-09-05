#include <fastcpd/fastcpd.h>

#include <armadillo>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;
using Row = std::unordered_map<std::string, std::string>;

std::vector<std::string> split(std::string const& line, char delimiter) {
  std::vector<std::string> fields;
  std::stringstream stream(line);
  std::string field;
  while (std::getline(stream, field, delimiter)) fields.push_back(field);
  if (!line.empty() && line.back() == delimiter) fields.emplace_back();
  return fields;
}

std::vector<Row> read_tsv(fs::path const& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("unable to open " + path.string());
  std::string line;
  if (!std::getline(input, line)) {
    throw std::runtime_error("empty TSV fixture: " + path.string());
  }
  std::vector<std::string> const header = split(line, '\t');
  std::vector<Row> rows;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::vector<std::string> const fields = split(line, '\t');
    if (fields.size() != header.size()) {
      throw std::runtime_error("invalid TSV row in " + path.string());
    }
    Row row;
    for (std::size_t index = 0; index < header.size(); ++index) {
      row.emplace(header[index], fields[index]);
    }
    rows.push_back(std::move(row));
  }
  return rows;
}

arma::mat read_csv(fs::path const& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("unable to open " + path.string());
  std::string line;
  if (!std::getline(input, line)) {
    throw std::runtime_error("empty CSV fixture: " + path.string());
  }
  std::size_t const columns = split(line, ',').size();
  std::vector<double> values;
  std::size_t rows = 0;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::vector<std::string> const fields = split(line, ',');
    if (fields.size() != columns) {
      throw std::runtime_error("invalid CSV row in " + path.string());
    }
    for (std::string const& field : fields) values.push_back(std::stod(field));
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

std::vector<double> parse_numbers(std::string const& value, char delimiter) {
  if (value.empty() || value == "-") return {};
  std::vector<double> result;
  for (std::string const& field : split(value, delimiter)) {
    if (field == "NA" || field == "NaN") {
      result.push_back(std::numeric_limits<double>::quiet_NaN());
    } else {
      result.push_back(std::stod(field));
    }
  }
  return result;
}

arma::colvec as_colvec(std::vector<double> const& values) {
  arma::colvec result(values.size());
  for (arma::uword index = 0; index < result.n_elem; ++index) {
    result(index) = values[index];
  }
  return result;
}

fastcpd::Options detector_options(Row const& row) {
  fastcpd::Options options;
  std::string const beta = row.at("beta");
  if (beta == "BIC" || beta == "MBIC" || beta == "MDL") {
    options.beta_criterion = beta;
  } else {
    options.beta = std::stod(beta);
  }
  options.cost_adjustment = row.at("cost_adjustment");
  options.trim = std::stod(row.at("trim"));
  options.vanilla_percentage = std::stod(row.at("vanilla_percentage"));
  options.order = as_colvec(parse_numbers(row.at("order"), ','));
  std::vector<double> const p_response =
      parse_numbers(row.at("p_response"), ',');
  if (!p_response.empty()) {
    options.p_response = static_cast<unsigned int>(p_response.front());
  }
  std::vector<double> const variance =
      parse_numbers(row.at("variance_estimation"), ',');
  if (!variance.empty()) {
    options.variance_estimate = arma::diagmat(as_colvec(variance));
  }
  std::vector<double> const seed =
      parse_numbers(row.at("random_state"), ',');
  if (!seed.empty()) options.seed = static_cast<std::int32_t>(seed.front());

  return options;
}

fastcpd::Result run_detector(Row const& row, arma::mat const& data) {
  fastcpd::Options options = detector_options(row);

  std::string const family = row.at("family");
  if (family == "mean") return fastcpd::detect_mean(data, options);
  if (family == "variance") return fastcpd::detect_variance(data, options);
  if (family == "meanvariance") {
    return fastcpd::detect_meanvariance(data, options);
  }
  if (family == "exponential") {
    return fastcpd::detect_exponential(data, options);
  }
  if (family == "lm") return fastcpd::detect_lm(data, options);
  if (family == "lasso") return fastcpd::detect_lasso(data, options);
  if (family == "binomial") return fastcpd::detect_binomial(data, options);
  if (family == "poisson") return fastcpd::detect_poisson(data, options);
  if (family == "quantile") return fastcpd::detect_quantile(data, options);
  if (family == "var") return fastcpd::detect_var(data, options);
  if (family == "rank") return fastcpd::detect_rank(data, options);
  if (family == "kcp") return fastcpd::detect_kernel(data, options);
  if (family == "ar") return fastcpd::detect_ar(data.col(0), options);
  if (family == "arma") return fastcpd::detect_arma(data.col(0), options);
  if (family == "arima") return fastcpd::detect_arima(data.col(0), options);
  if (family == "garch") return fastcpd::detect_garch(data.col(0), options);
  throw std::runtime_error("unregistered detector family: " + family);
}

struct ExpectedField {
  std::vector<std::size_t> shape;
  std::vector<double> values;
  double tolerance;
};

using ExpectedCase = std::unordered_map<std::string, ExpectedField>;

std::unordered_map<std::string, ExpectedCase> read_expected(
    fs::path const& path) {
  std::unordered_map<std::string, ExpectedCase> result;
  for (Row const& row : read_tsv(path)) {
    std::vector<double> const parsed_shape =
        parse_numbers(row.at("shape"), ',');
    ExpectedField field;
    for (double const value : parsed_shape) {
      field.shape.push_back(static_cast<std::size_t>(value));
    }
    field.values = parse_numbers(row.at("values"), ';');
    field.tolerance = std::stod(row.at("tolerance"));
    result[row.at("case_id")].emplace(row.at("field"), std::move(field));
  }
  return result;
}

void compare_value(std::string const& label, double actual, double expected,
                   double tolerance) {
  if (std::isnan(expected)) {
    if (!std::isnan(actual)) {
      throw std::runtime_error(label + " expected NaN");
    }
    return;
  }
  if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
    std::ostringstream message;
    message.precision(17);
    message << label << " differs: actual=" << actual
            << " expected=" << expected << " tolerance=" << tolerance;
    throw std::runtime_error(message.str());
  }
}

void compare_vector(std::string const& label, arma::colvec const& actual,
                    ExpectedField const& expected) {
  if (expected.shape.size() != 1 || actual.n_elem != expected.shape.front() ||
      actual.n_elem != expected.values.size()) {
    throw std::runtime_error(label + " shape differs");
  }
  for (arma::uword index = 0; index < actual.n_elem; ++index) {
    compare_value(label, actual(index), expected.values[index],
                  expected.tolerance);
  }
}

void compare_matrix(std::string const& label, arma::mat const& actual,
                    ExpectedField const& expected) {
  if (expected.shape.size() != 2 || actual.n_rows != expected.shape[0] ||
      actual.n_cols != expected.shape[1] ||
      actual.n_elem != expected.values.size()) {
    throw std::runtime_error(label + " shape differs");
  }
  std::size_t index = 0;
  for (arma::uword row = 0; row < actual.n_rows; ++row) {
    for (arma::uword column = 0; column < actual.n_cols; ++column) {
      compare_value(label, actual(row, column), expected.values[index++],
                    expected.tolerance);
    }
  }
}

void compare_detector(std::string const& case_id, fastcpd::Result const& result,
                      ExpectedCase const& expected) {
  compare_vector(case_id + " cp_set", result.change_points,
                 expected.at("cp_set"));
  compare_vector(case_id + " raw_cp_set", result.raw_change_points,
                 expected.at("raw_cp_set"));
  compare_vector(case_id + " cost_values", result.cost_values,
                 expected.at("cost_values"));
  compare_matrix(case_id + " residuals", result.residuals,
                 expected.at("residuals"));
  compare_matrix(case_id + " thetas", result.thetas, expected.at("thetas"));
}

double confidence_value(fastcpd::ConfidenceInterval const& row,
                        std::string const& field) {
  if (field == "index") return row.index;
  if (field == "segment") return row.segment;
  if (field == "parameter") return row.parameter;
  if (field == "estimate") return row.estimate;
  if (field == "lower") return row.lower;
  if (field == "upper") return row.upper;
  if (field == "detection_rate") return row.detection_rate;
  if (field == "profile_min") return row.profile_min;
  if (field == "cutoff") return row.cutoff;
  if (field == "se") return row.se;
  if (field == "level") return row.level;
  throw std::runtime_error("unregistered confidence field: " + field);
}

void compare_confidence(
    std::string const& case_id,
    std::vector<fastcpd::ConfidenceInterval> const& actual,
    ExpectedCase const& expected) {
  for (auto const& [field, values] : expected) {
    if (values.shape.size() != 1 || values.shape.front() != actual.size() ||
        values.values.size() != actual.size()) {
      throw std::runtime_error(case_id + " " + field + " shape differs");
    }
    for (std::size_t index = 0; index < actual.size(); ++index) {
      compare_value(case_id + " " + field,
                    confidence_value(actual[index], field),
                    values.values[index], values.tolerance);
    }
  }
}

fastcpd::ConfidenceOptions confidence_options(
    Row const& row, Row const& source) {
  fastcpd::ConfidenceOptions options;
  std::string const operation = row.at("operation");
  options.parm = operation == "confint_wald" ? "theta" : "cp";
  if (operation == "confint_profile") options.method = "profile";
  if (operation == "confint_wald") options.method = "wald";
  if (operation == "confint_bootstrap") options.method = "bootstrap";
  options.level = std::stod(row.at("level"));
  std::vector<double> const replicates = parse_numbers(row.at("B"), ',');
  if (!replicates.empty()) {
    options.bootstrap_replicates =
        static_cast<unsigned int>(replicates.front());
  }
  std::vector<double> const window = parse_numbers(row.at("window"), ',');
  if (!window.empty()) {
    options.window = static_cast<unsigned int>(window.front());
  }
  std::vector<double> const seed =
      parse_numbers(row.at("random_state"), ',');
  if (!seed.empty()) options.seed = static_cast<std::int32_t>(seed.front());
  options.detector_options = detector_options(source);
  return options;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    fs::path const fixtures = find_fixture_directory(argc, argv);
    std::unordered_map<std::string, ExpectedCase> const expected =
        read_expected(fixtures / "expected_outputs.tsv");
    std::vector<Row> const manifest = read_tsv(fixtures / "manifest.tsv");
    std::unordered_map<std::string, Row> source_rows;
    std::unordered_map<std::string, arma::mat> source_data;
    std::unordered_map<std::string, fastcpd::Result> source_results;
    std::size_t detector_count = 0;
    for (Row const& row : manifest) {
      if (row.at("operation") != "detect") continue;
      std::string const case_id = row.at("case_id");
      arma::mat const data = read_csv(fixtures / row.at("data_file"));
      fastcpd::Result const result = run_detector(row, data);
      if (result.family != row.at("family")) {
        throw std::runtime_error(case_id + " public family metadata differs");
      }
      compare_detector(case_id, result, expected.at(case_id));
      source_rows.emplace(case_id, row);
      source_data.emplace(case_id, data);
      source_results.emplace(case_id, result);
      ++detector_count;
    }
    if (detector_count != 19) {
      throw std::runtime_error("unexpected detector fixture count");
    }
    std::size_t confidence_count = 0;
    for (Row const& row : manifest) {
      if (row.at("operation").rfind("confint_", 0) != 0) continue;
      std::string const case_id = row.at("case_id");
      std::string const source_case = row.at("source_case");
      std::vector<fastcpd::ConfidenceInterval> const intervals =
          fastcpd::confint(source_results.at(source_case),
                           source_data.at(source_case),
                           confidence_options(row, source_rows.at(source_case)));
      compare_confidence(case_id, intervals, expected.at(case_id));
      ++confidence_count;
    }
    if (confidence_count != 13) {
      throw std::runtime_error("unexpected confidence fixture count");
    }
    std::cout << detector_count << " standalone detector and "
              << confidence_count << " confidence fixture contracts passed\n";
    return 0;
  } catch (std::exception const& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
