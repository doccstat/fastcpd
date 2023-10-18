## Updates since last CRAN release (0.7.2)

### fastcpd 0.8.4

*   Add cheatsheets.
*   Refactor code and utilize the `cost_function_wrapper`.
*   Optimize warm start.

### fastcpd 0.8.3

*   Add `fastcpd.ts` / `fastcpd_ts` for time series data.
*   Fix pre segmengation bug for `lasso`.
*   Fix bug related to `vanilla_percentage` parameter for `lasso`.
*   Add tests with invalid family for `fastcpd.ts`.
*   Remove the `cp_only = TRUE` default when the family is "custom".
*   Improved plotting for "ar" and "var" families.
*   Add test coverage for `cp_only = TRUE` and `fastcpd_ts`.
*   Increase test coverage.
*   Provide user selection when `ggplot2` is not installed.

### fastcpd 0.8.2

*   Add cheatsheets WIP.
*   Add smaller examples test for penalized linear regression.

### fastcpd 0.8.1

*   Add new "ar" family for autoregressive models.
*   Add new "var" family for vector autoregressive models.

### fastcpd 0.8.0

*   Deal with the following:

        Due to the excessive calls to `glmnet` between R and C++,
        it is better to use the R implementation of `fastcpd` for lasso.

*   Separate the use of internal C++ cost functions and user-defined R cost
    functions.
*   Add Codecov Icicle plot in README.
*   Remove `cost_optim` and `cost_update` from `RcppExports.R`.
*   Estimate the variance in the "gaussian" family dynamically.

### fastcpd 0.7.6

*   Move default cost functions definition inside the `fastcpd` definition.
*   Define constant unordered set to store the family sets.
*   Avoid using `length(formals(cost))` to check the number of arguments of
    `cost` function.
*   Introduce an internal family "vanilla".

### fastcpd 0.7.5

*   Add variance estimation example in linear regression.
*   Update reference page.
*   Add validation for `family`.
*   Add user selection when `ggplot2` is not installed.
*   Add AR(1) using `forecast` example in the tests.

### fastcpd 0.7.4

*   Update website UI.
*   Update `fastcpd` documentation.

### fastcpd 0.7.3

*   Allow multiple response variables in the `formula`.
*   Add fastcpd logo to README.
