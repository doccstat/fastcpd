## Updates since last CRAN release (1.0.0)

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
