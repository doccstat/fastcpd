# Changelog

## 1.3.0

This is the first coordinated R, Python, and standalone C++ source contract.

- R, Python, and standalone C++ sources now report version 1.3.0. Release CI
  rejects mismatched tags, wheel/sdist filenames, or archive metadata before
  publication, and the Python package is classified as production/stable.
- Scalar integer seeds reproduce R's KCP features and bootstrap sampling;
  NumPy generator objects remain a native Python extension.
- Shared fixtures cover every portable detector family and the common
  bootstrap, profile, and Wald confidence paths, including raw change points,
  costs, residuals, parameters, and interval diagnostics.
- Residuals use an observation-by-response matrix and autoregressive lag rows
  are padded with `NaN`.
- `detect_var()` accepts raw unlagged series only. Use advanced
  `family="mgaussian"` for a constructed response/predictor matrix.
- Multivariate-LM Wald intervals are explicitly unsupported, and complete
  binomial separation returns undefined (`NaN`) Wald uncertainty.
- Callable custom costs and `multiple_epochs` schedules remain R-only so the
  Python detector stays GIL-free.
- Public reference pages include executable examples. The package no longer
  ships `py.typed`; runtime annotations are not yet a complete static typing
  contract.

## 0.23.0

- Support CPython 3.11 through 3.14 with tested wheels for Linux, macOS, and
  Windows.
- Provide the independent Python API baseline used for the 1.3.0 parity work.
