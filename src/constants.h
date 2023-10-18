#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <algorithm>
#include <array>

using ::std::array;

template<typename C, typename T>
bool array_contains(C&& c, T e) {
    return std::find(std::begin(c), std::end(c), e) != std::end(c);
};

constexpr array CUSTOM_FAMILIES = {"custom", "vanilla"};
constexpr array FASTCPD_FAMILIES = {"gaussian", "binomial", "poisson", "lasso"};

#endif  // CONSTANTS_H_
