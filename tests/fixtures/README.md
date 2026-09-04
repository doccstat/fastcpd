# Shared fastcpd fixtures

This directory contains deterministic, language-neutral inputs used by the R
and Python test suites and by the standalone C++ mean-change example. Decimal
values are committed directly so every language evaluates the same data
without depending on an R or NumPy random-number stream.

`manifest.tsv` is the manually maintained, authoritative test index. It names
each operation, its input CSV, its public arguments, and its expected result.
Use `-` where a field does not apply. The R and Python suites parse the same
rows independently; the C++ example selects the manifest's mean detector row.

| column | meaning |
| --- | --- |
| `case_id` | Stable name used by the test suites. |
| `data_file` | CSV file relative to this directory. |
| `operation` | `detect`, `estimate_variance`, or `estimate_variance_arma`. |
| `family` | Public detector or variance-estimator family. |
| `order` | Scalar order or comma-separated model order. |
| `beta` | Numeric penalty or criterion name for detector rows. |
| `cost_adjustment` | Detector cost adjustment (`BIC`, `MBIC`, or `MDL`). |
| `trim` | Detector boundary/minimum-distance trim proportion. |
| `vanilla_percentage` | Fraction evaluated with exact PELT. |
| `expected_cp` | One-based change points separated by `;`. |
| `expected_value` | Expected scalar variance result. |
| `tolerance` | Numeric comparison tolerance for `expected_value`. |

`generate_fixtures.py` owns only the CSV contents. It records the deterministic
constructions: fixed step patterns, a repeating regression-error cycle,
integrated ARIMA increments, and an AR(1) recurrence for the ARMA variance
input. It also verifies that its CSV names match the files indexed by the
manifest. From the package root, run:

```sh
python tests/fixtures/generate_fixtures.py
python tests/fixtures/generate_fixtures.py --check
```

Do not hand-edit a generated CSV. Edit the construction in the generator,
regenerate the files, and update the manifest only when the operation or
expected result changes. Multiple manifest cases may intentionally reuse one
CSV. The Python source distribution includes the manifest, generator, README,
and all fixture CSVs through the fixture glob in `pyproject.toml`.

Response columns precede predictor columns, and change-point indices use the
one-based convention exposed by the public APIs. No language-specific
serializer or random seed is part of this contract.

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
