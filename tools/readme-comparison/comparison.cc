// Native README benchmark: run as readme_comparison observations.txt.
#include <fastcpd/fastcpd.h>
#include "colibri.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

int main(int argc, char** argv) {
  try {
    if (argc != 2) throw std::runtime_error("Supply the observation file");
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("Cannot open observation file");
    std::vector<double> x;
    for (double value; input >> value;) x.push_back(value);
    if (!input.eof() || x.size() != 1000000 ||
        !std::all_of(x.begin(), x.end(), [](double v) { return std::isfinite(v); })) {
      throw std::runtime_error("Expected 1000000 finite observations");
    }
    const int n = static_cast<int>(x.size());
    const arma::mat data(x.data(), x.size(), 1);
    auto run_fastcpd = [&]() {
      fastcpd::Options options;
      options.beta_criterion = "MBIC";
      options.cost_adjustment = "MBIC";
      options.variance_estimate = arma::eye<arma::mat>(1, 1);
      options.cp_only = true;
      auto result = fastcpd::detect_mean(data, options);
      return std::vector<int>(result.change_points.begin(), result.change_points.end());
    };
    auto run_fpop = [&]() {
      const double penalty = 2 * std::log(n);
      const auto bounds = std::minmax_element(x.begin(), x.end());
      std::vector<int> previous(n);
      std::vector<double> costs(n);
      colibri_op_c(x.data(), &n, &penalty, &*bounds.first, &*bounds.second,
                   previous.data(), costs.data());
      std::vector<int> points;
      for (int end = previous[n - 1]; end > 0; end = previous[end - 1]) {
        points.push_back(end);
      }
      std::reverse(points.begin(), points.end());
      return points;
    };
    const std::vector<int> expected{n / 2};
    if (run_fastcpd() != expected || run_fpop() != expected) {
      throw std::runtime_error("Warm-up change points differ from expected 500000");
    }
    std::cout << std::fixed << std::setprecision(9);
    for (int repeat = 0; repeat < 10; ++repeat) {
      for (int turn = 0; turn < 2; ++turn) {
        const bool use_fastcpd = (repeat + turn) % 2 == 0;
        const auto start = std::chrono::steady_clock::now();
        const auto points = use_fastcpd ? run_fastcpd() : run_fpop();
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        if (points != expected) throw std::runtime_error("Unexpected change points");
        std::cout << (use_fastcpd ? "fastcpd" : "fpop") << '\t'
                  << points.front() << '\t' << seconds << '\n';
      }
    }
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
