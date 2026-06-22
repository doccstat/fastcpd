## Updates since last CRAN release (1.0.0)

### fastcpd 1.0.5

*   Add favicon.

### fastcpd 1.0.4

*   Update macOS CI runner from deprecated `macos-13` to `macos-latest`.

### fastcpd 1.0.3

*   Remove `glmnet` and `Matrix` from `Imports`; lasso fitting now uses a
    pure C++ coordinate descent implementation with 5-fold cross-validation
    for lambda selection, eliminating two heavy runtime dependencies.
*   Replace `Matrix::nearPD()` with a pure-R eigendecomposition-based
    nearest positive-definite projection (`nearest_pd_`).

### fastcpd 1.0.2

*   Use a fairer comparison in README.

### fastcpd 1.0.1

*   Add `exponential` family to Python package (`fastcpd.segmentation.exponential`).
*   Fix Bazel build for Python CI: add missing `fastcpd_optim.h` and
    `families/exponential.h` to `hdrs` in `src/BUILD.bazel`.
*   Improve README benchmark methodology: use subprocess isolation via
    `callr::r()` inside `microbenchmark` with baseline subtraction for fair
    large-scale (n = 10^8) algorithm comparison.
*   Skip slow tests during coverage runs to prevent CI timeout.
