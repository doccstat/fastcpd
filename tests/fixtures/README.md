# Shared fastcpd fixtures

This directory contains deterministic, language-neutral inputs used by the R,
Python, and standalone C++ contract suites. Decimal
values are committed directly so every language evaluates the same data
without depending on an R or NumPy random-number stream.

`manifest.tsv` is the manually maintained, authoritative test index. It names
each operation, its input CSV, its public arguments, and its expected result.
Use `-` where a field does not apply. Detector rows cover every portable
built-in family, and confidence rows refer to their fitted detector through
`source_case`. The R and Python suites parse the same rows independently. The
C++ suite dispatches every detector and confidence row, compares every
detailed detector result and numeric interval diagnostic, and uses its
interface contract to check every variance-estimator row.

| column | meaning |
| --- | --- |
| `case_id` | Stable name used by the test suites. |
| `data_file` | CSV file relative to this directory. |
| `operation` | `detect`, `confint_bootstrap`, `confint_profile`, `confint_wald`, `estimate_variance`, or `estimate_variance_arma`. |
| `source_case` | Detector case used by a confidence operation. |
| `family` | Public detector or variance-estimator family. |
| `order` | Scalar order or comma-separated model order. |
| `beta` | Numeric penalty or criterion name for detector rows. |
| `cost_adjustment` | Detector cost adjustment (`BIC`, `MBIC`, or `MDL`). |
| `trim` | Detector boundary/minimum-distance trim proportion. |
| `vanilla_percentage` | Fraction evaluated with exact PELT. |
| `p_response` | Response-column count for multivariate linear regression. |
| `variance_estimation` | Optional scalar or comma-separated covariance diagonal. |
| `random_state` | Shared scalar seed for KCP or bootstrap. |
| `level`, `B`, `window` | Confidence-method controls. |
| `expected_cp` | One-based change points separated by `;`. |
| `expected_value` | Expected scalar variance result. |
| `tolerance` | Absolute numeric comparison tolerance for `expected_value`. |

`expected_outputs.tsv` is the normalized numerical contract generated from
the R reference implementation. Each row stores one vector or matrix field in
row-major order, including its exact shape and absolute comparison tolerance.
Detector rows cover `cp_set`, `raw_cp_set`, costs, residuals, and parameters;
confidence rows cover every numeric interval column. Multivariate residuals
are compared as conceptual `(observation, response)` matrices, independent of
R's current flattened S4 storage and Python's native two-dimensional storage.

KCP random features can put one bootstrap refit exactly at a native-math
decision boundary. The `kcp_bootstrap` detection rate therefore permits one
of its `B` detection indicators to differ across operating systems, plus
floating-point roundoff in the ratio. Its change point, interval endpoints,
and every other fixture field remain subject to their original tolerances.

`generate_fixtures.py` owns only the CSV contents. It records the deterministic
constructions: fixed step patterns, a repeating regression-error cycle,
integrated ARIMA increments, and an AR(1) recurrence for the ARMA variance
input. It also verifies that its CSV names match the files indexed by the
manifest. From the package root, run:

```sh
python tests/fixtures/generate_fixtures.py
python tests/fixtures/generate_fixtures.py --check
Rscript tests/fixtures/generate_r_shared_outputs.R
Rscript tests/fixtures/generate_r_shared_outputs.R --check
```

Do not hand-edit a generated CSV. Edit the construction in the generator,
regenerate the files, and update the manifest only when the operation or
expected result changes. Multiple manifest cases may intentionally reuse one
CSV. The Python source distribution includes the manifest, generator, README,
and all fixture CSVs through the fixture glob in `pyproject.toml`.

Response columns precede predictor columns, and change-point indices use the
one-based convention exposed by the public APIs. Scalar KCP/bootstrap seeds
use the explicitly documented R-compatible stream; native NumPy generator
objects remain a Python-specific extension and are outside these fixtures.

`r_rng.tsv` is a separate stochastic-conformance table generated from base R
with its default Mersenne-Twister, inversion-normal, and rejection-sampling
settings. It freezes the small subset of `set.seed()` behavior needed for
same-seed KCP features and confidence-bootstrap resampling. It is not indexed
by `manifest.tsv` because it describes an RNG stream rather than a detector
input.
The accompanying seeded feature matrices live under `stochastic/` so the
detector fixture generator's root-level CSV inventory remains exact.

Regenerate and verify the stochastic tables and seeded KCP feature matrix with:

```sh
Rscript tests/fixtures/generate_r_stochastic_fixtures.R
Rscript tests/fixtures/generate_r_stochastic_fixtures.R --check
```

The shared-output `--check` mode compares the current R results with the
committed rows using their absolute tolerances, so accepted BLAS/libm
roundoff does not rewrite the reference values. Running without `--check`
deliberately regenerates those values from the current R environment. An
output row may use a narrowly scoped tolerance above its manifest case's
default when only that field needs a cross-toolchain allowance.
