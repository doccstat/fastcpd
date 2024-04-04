#ifndef FASTCPD_CONSTANTS_H_
#define FASTCPD_CONSTANTS_H_

#include <algorithm>
#include <array>

#if defined(__cplusplus) && __cplusplus >= 201703L

template<typename C, typename T>
bool contain(C&& c, T e) {
  return std::find(std::begin(c), std::end(c), e) != std::end(c);
};
constexpr std::array FASTCPD_FAMILIES = {
  "gaussian", "binomial", "poisson", "lasso", "mgaussian",
  "arma",
  "mean", "variance", "meanvariance"
};

#else

#include <string>
#include <vector>
inline bool contain(std::vector<std::string> c, std::string e) {
  return std::find(std::begin(c), std::end(c), e) != std::end(c);
};
const std::vector<std::string> FASTCPD_FAMILIES = {
  "gaussian", "binomial", "poisson", "lasso", "mgaussian",
  "arma",
  "mean", "variance", "meanvariance"
};

#endif  // defined(__cplusplus) && __cplusplus >= 201703L

#endif  // FASTCPD_CONSTANTS_H_
