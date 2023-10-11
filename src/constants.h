#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <string>
#include <unordered_set>

using ::std::string;
using ::std::unordered_set;

// TODO(doccstat): `unordered_set` is not literal type and thus `constexpr`
// can not be used here. Blocked by `abseil` due to the non-header only issue.
const unordered_set<string> CUSTOM_FAMILIES = {"custom", "vanilla"};

const unordered_set<string> FASTCPD_FAMILIES =
  {"gaussian", "binomial", "poisson", "lasso"};

#endif  // CONSTANTS_H_
