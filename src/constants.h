#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <algorithm>
#include <array>

#if defined(__cplusplus) && __cplusplus >= 201703L
template<typename C, typename T>
bool contain(C&& c, T e) {
  return std::find(std::begin(c), std::end(c), e) != std::end(c);
};
#else
bool contain(std::array, std::string) {
  return std::find(std::begin(c), std::end(c), e) != std::end(c);
};
#endif

constexpr std::array CUSTOM_FAMILIES = {"custom", "vanilla"};
constexpr std::array FASTCPD_FAMILIES = {"gaussian", "binomial", "poisson", "lasso"};

#endif  // CONSTANTS_H_
