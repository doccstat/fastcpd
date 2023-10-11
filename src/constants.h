#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <string>
#include <unordered_set>

// TODO(doccstat): `std::unordered_set` is not literal type and thus `constexpr`
// can not be used here. Blocked by `abseil` due to the non-header only issue.
const std::unordered_set<std::string> CUSTOM_FAMILIES = {
  "custom", "vanilla"
};

const std::unordered_set<std::string> FASTCPD_FAMILIES = {
  "gaussian", "binomial", "poisson", "lasso",
};

#endif  // CONSTANTS_H
