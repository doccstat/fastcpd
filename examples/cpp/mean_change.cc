#include <fastcpd/fastcpd.h>

#include <armadillo>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

namespace fs = std::filesystem;

std::vector<std::string> split(std::string const& line, char delimiter) {
  std::vector<std::string> fields;
  std::stringstream stream(line);
  std::string field;
  while (std::getline(stream, field, delimiter)) fields.push_back(field);
  if (!line.empty() && line.back() == delimiter) fields.emplace_back();
  return fields;
}

struct FixtureCase {
  fs::path data_file;
  std::string family;
  double beta;
  std::string cost_adjustment;
  double trim;
  double vanilla_percentage;
  std::vector<double> expected_cp;
};

fs::path find_manifest(int argc, char** argv) {
  std::vector<fs::path> candidates;
  if (argc > 1) {
    fs::path supplied(argv[1]);
    candidates.push_back(supplied.filename() == "manifest.tsv"
                             ? supplied
                             : supplied / "manifest.tsv");
  } else {
    candidates = {
        fs::path("tests/fixtures/manifest.tsv"),
        fs::path("../tests/fixtures/manifest.tsv"),
        fs::path("../../tests/fixtures/manifest.tsv"),
    };
  }
  for (fs::path const& candidate : candidates) {
    std::error_code error;
    if (fs::is_regular_file(candidate, error)) return candidate;
  }
  throw std::runtime_error(
      "unable to locate tests/fixtures/manifest.tsv; pass its directory "
      "as the first argument");
}

FixtureCase read_mean_fixture_case(fs::path const& manifest_path) {
  std::ifstream input(manifest_path);
  if (!input) {
    throw std::runtime_error("unable to open fixture manifest: " +
                             manifest_path.string());
  }

  std::string line;
  if (!std::getline(input, line)) {
    throw std::runtime_error("fixture manifest is empty");
  }
  std::vector<std::string> header = split(line, '\t');
  std::unordered_map<std::string, std::size_t> column;
  for (std::size_t i = 0; i < header.size(); ++i) column[header[i]] = i;
  for (char const* required : {"data_file", "operation", "family", "beta",
                               "cost_adjustment", "trim",
                               "vanilla_percentage", "expected_cp"}) {
    if (column.find(required) == column.end()) {
      throw std::runtime_error("fixture manifest lacks column: " +
                               std::string(required));
    }
  }

  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::vector<std::string> fields = split(line, '\t');
    if (fields.size() < header.size()) continue;
    auto value = [&fields, &column](char const* name) -> std::string const& {
      return fields[column.at(name)];
    };
    if (value("operation") != "detect" || value("family") != "mean") {
      continue;
    }
    std::vector<double> expected;
    for (std::string const& item : split(value("expected_cp"), ';')) {
      if (!item.empty()) expected.push_back(std::stod(item));
    }
    return FixtureCase{
        manifest_path.parent_path() / value("data_file"),
        value("family"),
        std::stod(value("beta")),
        value("cost_adjustment"),
        std::stod(value("trim")),
        std::stod(value("vanilla_percentage")),
        std::move(expected),
    };
  }
  throw std::runtime_error("mean detector fixture not found");
}

arma::mat read_csv_column(fs::path const& data_path) {
  std::ifstream input(data_path);
  if (!input) {
    throw std::runtime_error("unable to open fixture data: " +
                             data_path.string());
  }
  std::string line;
  if (!std::getline(input, line)) {
    throw std::runtime_error("fixture data is empty: " + data_path.string());
  }
  std::vector<double> values;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::size_t separator = line.find(',');
    if (separator != std::string::npos) line.resize(separator);
    values.push_back(std::stod(line));
  }
  if (values.empty()) {
    throw std::runtime_error("fixture data has no observations: " +
                             data_path.string());
  }
  arma::mat data(values.size(), 1);
  for (arma::uword i = 0; i < data.n_rows; ++i) data(i, 0) = values[i];
  return data;
}

}  // namespace

int main(int argc, char** argv) {
  try {
    fs::path manifest = find_manifest(argc, argv);
    FixtureCase fixture = read_mean_fixture_case(manifest);
    arma::mat data = read_csv_column(fixture.data_file);

    fastcpd::Options options;
    options.family = fixture.family;
    options.beta = fixture.beta;
    options.cost_adjustment = fixture.cost_adjustment;
    options.trim = fixture.trim;
    options.vanilla_percentage = fixture.vanilla_percentage;
    options.cp_only = true;

    fastcpd::Result result = fastcpd::detect(data, options);

    result.change_points.t().print(std::cout, "change points:");
    bool matches = result.change_points.n_elem == fixture.expected_cp.size();
    for (arma::uword i = 0; matches && i < result.change_points.n_elem; ++i) {
      matches = std::abs(result.change_points(i) - fixture.expected_cp[i]) <=
                1e-9;
    }
    if (!matches) {
      std::cerr << "change points do not match " << manifest << "\n";
      return 1;
    }
    return 0;
  } catch (std::exception const& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
