## Updates since last CRAN release (1.0.0)

### fastcpd 1.1.0

*   Add `fastcpd_kcp()` / `fastcpd.kcp()`: kernel change-point detection via
    random Fourier feature approximation of the RBF kernel.
*   Add `fastcpd_rank()` / `fastcpd.rank()`: distribution-free rank-based
    change-point detection.
*   Replace `RProgress` dependency with a built-in progress bar; remove
    `progress` from `LinkingTo`. Rename `r.progress` → `show.progress`;
    default changed to `FALSE`.
*   Unify ARMA and MA C++ implementation: remove `#ifndef NO_RCPP` split;
    both families now use a native C++ Kalman ML path in R and Python.
*   Python: add `garch`, `lm`, `binomial`, `poisson`, `arma`, `ar`, `arima`
    to `fastcpd.segmentation` (all 13 non-custom families now supported).
*   Fix `NEWS.md` history: duplicate `fastcpd 0.99.9` entry corrected to
    `fastcpd 0.9.9`.

### fastcpd 1.0.8

*   Add quantile regression family (`fastcpd_quantile`): robust change-point
    detection via pinball loss with IRLS fitting; resistant to outliers unlike
    least-squares methods.

### fastcpd 1.0.7

*   Change default `trim` from `0.05` to `0`: change points near segment
    boundaries are no longer suppressed by default.
*   Reduce GARCH example sizes: `n` reduced from 1501/200/300 to 401/120/150
    for faster documentation builds and test runs.
*   Replace boxplot with violin plot in README benchmark figure.

### fastcpd 1.0.6

*   Performance improvements: absl::InlinedVector for pruned set (O(1) memory
    vs O(n) pre-allocation), cross-platform prefetch via absl::PrefetchToLocalCache,
    and baseline-subtracted benchmark plot with horizontal layout.

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
