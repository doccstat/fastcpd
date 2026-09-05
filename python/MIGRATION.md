# Migrating to the coordinated 1.3.0 interface

Python moves directly from the independent 0.x version line to 1.3.0 because
this release is the first source contract shared with R and standalone C++.
No placeholder Python releases are required for the skipped version numbers.

## VAR input

`detect_var(data, order=p)` now always interprets `data` as a raw unlagged
multivariate series. The old shape-based `p_response` heuristic was ambiguous
and is removed.

For a matrix whose first `q` columns are responses and remaining columns are
already-constructed predictors, use the advanced interface explicitly:

```python
result = fastcpd.detect(
    data=design,
    family="mgaussian",
    p_response=q,
    beta=penalty,
)
```

## Results and confidence intervals

- `raw_cp_set` is the untrimmed native traceback; `cp_set` remains the public
  trimmed/merged result.
- Multivariate residuals have shape `(n_observations, n_responses)`, with one
  leading all-`NaN` row per autoregressive lag.
- Multivariate linear-model Wald intervals now raise `NotImplementedError`
  until R and Python share a validated covariance estimator.
- Complete binomial separation returns `NaN` standard errors and bounds.
- Rank fits retain original observations in Python's `result.data`; fitted
  costs, residuals, parameters, and intervals still use centered ranks.

## Randomness and callbacks

Scalar integer `random_state` values use the R-compatible stream for KCP and
bootstrap parity. Passing a NumPy `Generator` or `RandomState` keeps native
NumPy semantics.

Python continues to reject callable custom costs and `multiple_epochs`
schedules. R and standalone C++ expose language-native callback extensions;
omitting a Python callback ABI preserves the GIL-free native detector path.

## Static typing

The previous empty `py.typed` marker overstated the package's annotation
coverage. It has been removed for 1.3.0. Runtime signatures remain stable, but
consumers should not treat the package as a complete PEP 561 typed library.
