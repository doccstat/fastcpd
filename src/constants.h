#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <algorithm>

#if defined(__cplusplus) && __cplusplus >= 201703L

#include <array>
template<typename C, typename T>
bool contain(C&& c, T e) {
  return std::find(std::begin(c), std::end(c), e) != std::end(c);
};
constexpr std::array CUSTOM_FAMILIES = {"custom", "vanilla"};
constexpr std::array FASTCPD_FAMILIES =
  {"gaussian", "binomial", "poisson", "lasso"};

#else

#include <string>
#include <vector>
inline bool contain(std::vector<std::string> c, std::string e) {
  return std::find(std::begin(c), std::end(c), e) != std::end(c);
};
const std::vector<std::string> CUSTOM_FAMILIES = {"custom", "vanilla"};
const std::vector<std::string> FASTCPD_FAMILIES =
  {"gaussian", "binomial", "poisson", "lasso"};

#endif  // defined(__cplusplus) && __cplusplus >= 201703L

#endif  // CONSTANTS_H_
