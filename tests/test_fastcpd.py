import csv
import inspect
import subprocess
import sys
import threading
import time
import unittest
from pathlib import Path

import fastcpd as fastcpd_pkg
import fastcpd.segmentation as segmentation
import numpy as np
import pytest
from fastcpd.interface import fastcpd_impl as native_fastcpd_impl
from fastcpd.segmentation import (
    ar, arima, arma, binomial, detect, detect_kernel, detect_mean,
    detect_quantile, detect_rank, exponential, garch, lasso, lm, mean,
    meanvariance, poisson, var, variance,
)
from numpy import concatenate
from numpy.random import exponential as rexp
from numpy.random import multivariate_normal, randn, seed


_SHARED_FIXTURE_DIR = Path(__file__).resolve().parent / 'fixtures'
_SHARED_MANIFEST_COLUMNS = (
    'case_id', 'data_file', 'operation', 'family', 'order', 'beta',
    'cost_adjustment', 'trim', 'vanilla_percentage', 'expected_cp',
    'expected_value', 'tolerance',
)


def _shared_fixture_manifest():
    """Read the language-neutral fixture manifest.

    The R suite reads the same TSV independently.  Keeping this loader in the
    test module (rather than in the package) ensures fixture files never become
    part of the runtime API or wheel payload.
    """
    with (_SHARED_FIXTURE_DIR / 'manifest.tsv').open(
        newline='', encoding='utf-8'
    ) as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        if tuple(reader.fieldnames or ()) != _SHARED_MANIFEST_COLUMNS:
            raise AssertionError('Unexpected shared fixture manifest columns')
        rows = list(reader)
    case_ids = [row['case_id'] for row in rows]
    if len(case_ids) != len(set(case_ids)):
        raise AssertionError('Shared fixture case_id values must be unique')
    return rows


def _shared_fixture_data(filename):
    """Load one shared fixture as a two-dimensional floating-point matrix."""
    return np.loadtxt(
        _SHARED_FIXTURE_DIR / filename,
        delimiter=',',
        skiprows=1,
        ndmin=2,
    )


def _shared_fixture_order(value):
    """Parse the manifest's scalar/comma-separated order representation."""
    values = tuple(int(item) for item in value.split(',') if item)
    return values[0] if len(values) == 1 else values


def _run_shared_fixture(row, data):
    """Run a manifest row through the corresponding public Python wrapper."""
    assert row['operation'] == 'detect'
    family = row['family']
    order = _shared_fixture_order(row['order'])
    beta = row['beta']
    if beta not in {'BIC', 'MBIC', 'MDL'}:
        beta = float(beta)
    kwargs = {
        'beta': beta,
        'cost_adjustment': row['cost_adjustment'],
        'trim': float(row['trim']),
        'vanilla_percentage': float(row['vanilla_percentage']),
    }
    if family == 'mean':
        return detect_mean(data, **kwargs)
    if family == 'variance':
        return variance(data, **kwargs)
    if family == 'meanvariance':
        return meanvariance(data, **kwargs)
    if family == 'exponential':
        return exponential(data, **kwargs)
    if family == 'lm':
        return lm(data, **kwargs)
    if family == 'var':
        return var(data, order=order, **kwargs)
    if family == 'rank':
        return detect_rank(data, **kwargs)
    if family == 'arima':
        return arima(data, order=order, **kwargs)
    raise AssertionError(f'No Python wrapper registered for {family!r}')


def _expit(x):
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    positive = x >= 0
    out[positive] = 1.0 / (1.0 + np.exp(-x[positive]))
    exp_x = np.exp(x[~positive])
    out[~positive] = exp_x / (1.0 + exp_x)
    return out


class TestBasic(unittest.TestCase):

    def test_shared_cpp_mean_contract(self):
        data = concatenate((np.zeros(50), np.full(50, 5.0)))
        result = detect_mean(
            data,
            beta=5.0,
            cost_adjustment='BIC',
            trim=0.0,
            variance_estimation=np.eye(1),
            cp_only=True,
        )
        self.assertIsInstance(result, segmentation.CpdResult)
        np.testing.assert_array_equal(result.cp_set, [50])
        self.assertEqual(result.cp_set.dtype, np.dtype(np.int64))
        self.assertFalse(result.details_available)
        self.assertEqual(result.cost_values.shape, (0,))
        self.assertEqual(result.residuals.shape, (0, 0))
        self.assertEqual(result.thetas.shape, (0, 0))
        self.assertEqual(result.data.shape, (100, 1))
        self.assertEqual(result.family, 'mean')
        self.assertEqual(result.order, (0, 0, 0))
        self.assertFalse(result.cp_set.flags.writeable)
        self.assertFalse(result.data.flags.writeable)

    def test_native_numpy_binding_releases_gil(self):
        data = np.concatenate([np.zeros(5000), np.ones(5000)]).reshape(-1, 1)
        worker_ran = threading.Event()

        def worker():
            time.sleep(0.05)
            worker_ran.set()

        thread = threading.Thread(target=worker)
        thread.start()
        native_result = native_fastcpd_impl(
            5.0, 'BIC', True, data, 1e-10, 'mean', np.array([1.0]),
            np.empty(0), 0.0, np.array([0.0, 0.0, 0.0]), 0, 0,
            float('-inf'), 2, 0.0, np.empty(0), 1.0, np.eye(1), False,
            False,
        )
        worker_ran_during_detection = worker_ran.is_set()
        thread.join()

        self.assertTrue(worker_ran_during_detection)
        self.assertIsInstance(native_result['cp_set'], np.ndarray)
        self.assertEqual(native_result['cp_set'].dtype, np.dtype(np.int64))
        self.assertIsInstance(native_result['residuals'], np.ndarray)

    def test_unified_interface_aliases(self):
        self.assertIs(fastcpd_pkg.detect_mean, fastcpd_pkg.mean)
        self.assertIs(fastcpd_pkg.detect_kernel, fastcpd_pkg.kernel)
        self.assertIs(fastcpd_pkg.detect_kcp, fastcpd_pkg.kcp)
        self.assertIs(fastcpd_pkg.detect_kcp, fastcpd_pkg.detect_kernel)
        self.assertIs(fastcpd_pkg.detect_lm, fastcpd_pkg.lm)
        self.assertIs(
            fastcpd_pkg.detect_mean_variance,
            fastcpd_pkg.detect_meanvariance,
        )
        self.assertIs(
            fastcpd_pkg.detect_linear_regression,
            fastcpd_pkg.detect_lm,
        )
        self.assertIs(
            fastcpd_pkg.detect_logistic_regression,
            fastcpd_pkg.detect_binomial,
        )
        self.assertIs(
            fastcpd_pkg.detect_poisson_regression,
            fastcpd_pkg.detect_poisson,
        )
        self.assertIs(
            fastcpd_pkg.detect_quantile_regression,
            fastcpd_pkg.detect_quantile,
        )
        self.assertIs(
            fastcpd_pkg.estimate_variance_mean,
            fastcpd_pkg.variance_estimation.mean,
        )
        for name in (
            'detect_time_series', 'detect_ts', 'fastcpd_ts', 'time_series', 'ts'
        ):
            self.assertFalse(hasattr(fastcpd_pkg, name))
            self.assertFalse(hasattr(segmentation, name))

        seed(17)
        data = concatenate((np.random.normal(0, 0.2, 40),
                            np.random.normal(3, 0.2, 40)))
        np.testing.assert_array_equal(
            detect_mean(data).cp_set, mean(data).cp_set
        )
        np.testing.assert_array_equal(
            detect_rank(data).cp_set,
            fastcpd_pkg.rank(data).cp_set,
        )
        np.testing.assert_array_equal(
            detect_kernel(
                data, order=(20, 1), random_state=17,
            ).cp_set,
            fastcpd_pkg.kcp(
                data, order=(20, 1), random_state=17,
            ).cp_set,
        )

        x = np.arange(80, dtype=float)
        regression_data = np.column_stack([
            concatenate((x[:40], -x[40:])),
            x,
        ])
        np.testing.assert_array_equal(
            fastcpd_pkg.detect_linear_regression(
                regression_data, beta=5.0
            ).cp_set,
            fastcpd_pkg.lm(regression_data, beta=5.0).cp_set,
        )

    def test_unified_variance_interface(self):
        data = np.array([0.0, 1.0, 2.0, 4.0])
        expected = fastcpd_pkg.variance_estimation.mean(data)
        np.testing.assert_allclose(
            fastcpd_pkg.estimate_variance_mean(data),
            expected,
        )
        np.testing.assert_allclose(
            fastcpd_pkg.estimate_variance(data, family='mean'),
            expected,
        )

    def test_detect_defaults_match_r(self):
        parameters = inspect.signature(detect).parameters
        self.assertEqual(parameters['trim'].default, 0.0)
        self.assertEqual(parameters['vanilla_percentage'].default, 0.0)
        self.assertIsNone(parameters['multiple_epochs'].default)

    def test_detect_rejects_unsupported_arguments(self):
        with self.assertRaisesRegex(TypeError, "unexpected keyword"):
            detect_mean(np.arange(10), unknown_option=True)
        with self.assertRaisesRegex(NotImplementedError, "multiple_epochs"):
            detect_mean(np.arange(10), multiple_epochs=lambda _: 1)
        with self.assertRaisesRegex(NotImplementedError, "intentionally R-only"):
            detect(
                data=np.arange(10), family='custom',
                cost=lambda segment: float(np.sum(segment)),
            )

    def test_quantile_interface(self):
        seed(18)
        x = randn(120)
        y = concatenate((2 * x[:60], -2 * x[60:])) + 0.05 * randn(120)
        result = detect_quantile(
            np.column_stack([y, x]), order=0.5, trim=0.05,
            vanilla_percentage=1.0,
        )
        self.assertGreater(len(result.cp_set), 0)
        self.assertAlmostEqual(result.cp_set[0], 60, delta=10)
        interval = result.confint(method='profile', level=0.8, window=2)
        self.assertEqual(interval[0]['estimate'], result.cp_set[0])
        self.assertLessEqual(interval[0]['lower'], result.cp_set[0])
        self.assertGreaterEqual(interval[0]['upper'], result.cp_set[0])

    def test_quantile_order_sequence_spellings(self):
        """One-element sequence orders match R's wrapper normalization."""
        rng = np.random.default_rng(180)
        data = np.column_stack([
            rng.normal(size=40), rng.normal(size=40)
        ])
        for order in ((0.5,), [0.5], np.array([0.5])):
            result = detect_quantile(
                data, order=order, beta=1e6, cp_only=True,
                vanilla_percentage=1.0,
            )
            self.assertEqual(result.order, (0.5,))
        with self.assertRaisesRegex(ValueError, 'one level'):
            detect_quantile(data, order=(0.25, 0.75), cp_only=True)

    def test_mean(self):
        seed(0)
        covariance_mat = [[100, 0, 0], [0, 100, 0], [0, 0, 100]]
        data = concatenate((multivariate_normal([0, 0, 0], covariance_mat, 300),
                            multivariate_normal(
                                [50, 50, 50], covariance_mat, 400),
                            multivariate_normal([2, 2, 2], covariance_mat, 300)
                            ))
        result = mean(data)
        self.assertEqual(result.cp_set[0], 300)
        self.assertEqual(result.cp_set[1], 700)

    def test_mean_confint(self):
        seed(16)
        data = concatenate((np.random.normal(0, 0.2, 40),
                            np.random.normal(3, 0.2, 40)))
        result = mean(data)
        interval = result.confint(
            method='profile', level=0.8, window=8
        )
        self.assertEqual(interval[0]['estimate'], result.cp_set[0])
        self.assertLessEqual(interval[0]['lower'], result.cp_set[0])
        self.assertGreaterEqual(interval[0]['upper'], result.cp_set[0])
        bootstrap_interval = result.confint(
            method='bootstrap', level=0.8, B=3, random_state=16,
            detect_kwargs={'beta': 5.0, 'cost_adjustment': 'BIC'},
        )
        self.assertEqual(bootstrap_interval[0]['estimate'], result.cp_set[0])

    def test_exponential(self):
        seed(1)
        data = concatenate((rexp(scale=1.0, size=500), rexp(scale=5.0, size=500)))
        result = exponential(data)
        self.assertEqual(result.cp_set[0], 504)

    def test_variance(self):
        seed(2)
        data = concatenate((np.random.normal(0, 1, 500), np.random.normal(0, 5, 500)))
        result = variance(data)
        self.assertEqual(result.cp_set[0], 501)

    def test_meanvariance(self):
        seed(3)
        data = concatenate((np.random.normal(0, 1, 300), np.random.normal(5, 3, 300)))
        result = meanvariance(data, trim=0.05)
        self.assertEqual(result.cp_set[0], 300)

    def test_var_mgaussian(self):
        # VAR(1) wrapper accepts the raw two-column time series.
        seed(4)
        q = 2
        cov = [[1, 0], [0, 1]]
        y_raw = concatenate((
            multivariate_normal([0, 0], cov, 300),
            multivariate_normal([5, 5], cov, 300),
        ))
        result = var(y_raw, order=1, trim=0.05)

        # Advanced mgaussian calls still accept a pre-constructed design.
        data_mg = np.column_stack([y_raw[1:], y_raw[:-1]])
        direct = detect(
            data=data_mg, family='mgaussian', p_response=q, order=(1,),
            trim=0.05,
        )
        self.assertGreater(len(direct.cp_set), 0)
        np.testing.assert_array_equal(
            result.cp_set,
            [change_point + 1 for change_point in direct.cp_set],
        )
        legacy = var(data_mg, order=1, p_response=q, trim=0.05)
        np.testing.assert_array_equal(legacy.cp_set, direct.cp_set)

    def test_lasso(self):
        seed(7)
        n, p = 400, 5
        X = randn(n, p)
        y1 = X[:200] @ np.array([3.0, 0, 0, 0, 0]) + randn(200) * 0.1
        y2 = X[200:] @ np.array([0, 0, 0, 0, -3.0]) + randn(200) * 0.1
        data = np.column_stack([concatenate([y1, y2]), X])
        result = lasso(data)
        self.assertEqual(result.cp_set[0], 200)

    def test_lm(self):
        seed(8)
        n = 400
        X = randn(n, 3)
        y = concatenate([
            X[:200] @ np.array([1.0, 0.0, 0.0]) + randn(200) * 0.5,
            X[200:] @ np.array([0.0, 0.0, 1.0]) + randn(200) * 0.5,
        ])
        data = np.column_stack([y, X])
        result = lm(data)
        self.assertGreater(len(result.cp_set), 0)
        self.assertAlmostEqual(result.cp_set[0], 200, delta=10)

    def test_binomial(self):
        seed(9)
        n = 600
        X = randn(n, 2)
        p1 = _expit(X[:300] @ np.array([2.0, 0.0]))
        p2 = _expit(X[300:] @ np.array([0.0, 2.0]))
        y = concatenate([
            np.random.binomial(1, p1),
            np.random.binomial(1, p2),
        ]).astype(float)
        data = np.column_stack([y, X])
        result = binomial(data)
        self.assertGreater(len(result.cp_set), 0)
        self.assertAlmostEqual(result.cp_set[0], 300, delta=20)
        profile = result.confint(method='profile', level=0.8, window=1)
        self.assertEqual(profile[0]['estimate'], result.cp_set[0])
        wald = result.confint(parm='theta', method='wald', level=0.8)
        self.assertEqual(len(wald), result.thetas.size)
        self.assertTrue(all(np.isfinite(row['se']) for row in wald))

    def test_poisson(self):
        seed(10)
        n = 600
        X = randn(n, 2)
        mu1 = np.exp(X[:300] @ np.array([0.8, 0.0]))
        mu2 = np.exp(X[300:] @ np.array([0.0, 0.8]))
        y = concatenate([
            np.random.poisson(mu1),
            np.random.poisson(mu2),
        ]).astype(float)
        data = np.column_stack([y, X])
        result = poisson(data, trim=0.05, vanilla_percentage=1.0)
        self.assertGreater(len(result.cp_set), 0)
        self.assertAlmostEqual(result.cp_set[0], 300, delta=20)
        profile = result.confint(method='profile', level=0.8, window=1)
        self.assertEqual(profile[0]['estimate'], result.cp_set[0])
        wald = result.confint(parm='theta', method='wald', level=0.8)
        self.assertEqual(len(wald), result.thetas.size)
        self.assertTrue(all(np.isfinite(row['se']) for row in wald))

    @pytest.mark.long
    def test_garch(self):
        seed(11)
        n = 600
        # Two GARCH(1,1) segments with very different persistence.
        from numpy.random import default_rng
        rng = default_rng(11)
        x = np.zeros(n)
        h = np.ones(n)
        # Segment 1: low volatility  α=0.05, β=0.10
        for t in range(1, 300):
            h[t] = 0.5 + 0.05 * x[t - 1] ** 2 + 0.10 * h[t - 1]
            x[t] = rng.normal(0, np.sqrt(h[t]))
        # Segment 2: high volatility α=0.30, β=0.60
        for t in range(300, n):
            h[t] = 0.5 + 0.30 * x[t - 1] ** 2 + 0.60 * h[t - 1]
            x[t] = rng.normal(0, np.sqrt(h[t]))
        result = garch(x, order=(1, 1))
        self.assertGreater(len(result.cp_set), 0)
        self.assertAlmostEqual(result.cp_set[0], 300, delta=30)

    def test_ar(self):
        seed(12)
        n = 600
        x = np.zeros(n)
        # Segment 1: AR(1) with φ=0.8
        for t in range(1, 300):
            x[t] = 0.8 * x[t - 1] + randn()
        # Segment 2: AR(1) with φ=-0.8
        for t in range(300, n):
            x[t] = -0.8 * x[t - 1] + randn()
        result = ar(x, order=1)
        self.assertGreater(len(result.cp_set), 0)
        self.assertAlmostEqual(result.cp_set[0], 300, delta=20)

    def test_arma(self):
        seed(13)
        n = 600
        x = np.zeros(n)
        eps = randn(n)
        # Segment 1: ARMA(1,1) with φ=0.5, θ=0.3
        for t in range(1, 300):
            x[t] = 0.5 * x[t - 1] + eps[t] + 0.3 * eps[t - 1]
        # Segment 2: ARMA(1,1) with φ=-0.5, θ=-0.3
        eps2 = randn(n)
        for t in range(300, n):
            x[t] = -0.5 * x[t - 1] + eps2[t] - 0.3 * eps2[t - 1]
        result = arma(
            x, order=(1, 1), trim=0.05, vanilla_percentage=1.0,
        )
        self.assertGreater(len(result.cp_set), 0)
        self.assertAlmostEqual(result.cp_set[0], 300, delta=30)

    def test_arima(self):
        # Shared R/Python contract: difference each candidate segment, retain
        # original-series indices, and omit the cross-boundary difference.
        small = np.tile([0.1, -0.1], 20)
        large = np.resize(np.array([2.0, -2.0]), 41)
        x = np.concatenate([
            [0.0], np.cumsum(np.concatenate([small, large]))
        ])
        result = arima(x, order=(0, 1, 0))

        np.testing.assert_array_equal(result.cp_set, [41])
        mbic_adjustment = np.log(41) / 2
        expected_costs = [
            20 * (np.log(2 * np.pi) + np.log(0.01) + 1) + mbic_adjustment,
            20 * (np.log(2 * np.pi) + np.log(4.0) + 1) + mbic_adjustment,
        ]
        np.testing.assert_allclose(result.cost_values, expected_costs)
        np.testing.assert_allclose(result.thetas, [[0.01, 4.0]])

        residuals = np.asarray(result.residuals)[:, 0]
        np.testing.assert_array_equal(
            np.flatnonzero(np.isnan(residuals)) + 1, [1, 42]
        )
        np.testing.assert_allclose(
            residuals[~np.isnan(residuals)],
            np.concatenate([small, large[1:]]),
        )
        with self.assertRaisesRegex(ValueError, "include_mean=True"):
            arima(x, order=(0, 1, 0), include_mean=True)
        interval = result.confint(method='profile', level=0.8, window=1)
        self.assertEqual(interval[0]['estimate'], 41)
        self.assertLessEqual(interval[0]['lower'], 41)
        self.assertGreaterEqual(interval[0]['upper'], 41)

    def test_arima_d0_matches_arma(self):
        # ARIMA(1, 0, 0) should produce the same result as arma(x, order=(1, 0)).
        seed(15)
        n = 400
        x = np.zeros(n)
        for t in range(1, 200):
            x[t] = 0.8 * x[t - 1] + randn()
        for t in range(200, n):
            x[t] = -0.8 * x[t - 1] + randn()
        r1 = arima(x, order=(1, 0, 0))
        r2 = arma(x, order=(1, 0))
        np.testing.assert_array_equal(r1.cp_set, r2.cp_set)


_SHARED_ROWS = {row['case_id']: row for row in _shared_fixture_manifest()}
_SHARED_DETECTOR_ROWS = {
    case_id: row for case_id, row in _SHARED_ROWS.items()
    if row['operation'] == 'detect'
}
_SHARED_GENERIC_DETECTOR_ROWS = {
    case_id: row for case_id, row in _SHARED_DETECTOR_ROWS.items()
    if row['family'] != 'rank'
}
_SHARED_VARIANCE_ROWS = {
    case_id: row for case_id, row in _SHARED_ROWS.items()
    if row['operation'].startswith('estimate_variance')
}


@pytest.mark.parametrize(
    'case_id',
    tuple(_SHARED_DETECTOR_ROWS),
    ids=tuple(_SHARED_DETECTOR_ROWS),
)
def test_shared_fixture_change_point_contract(case_id):
    """The deterministic fixture manifest defines a cross-language CP contract."""
    row = _SHARED_DETECTOR_ROWS[case_id]
    data = _shared_fixture_data(row['data_file'])
    result = _run_shared_fixture(row, data)
    expected = np.asarray(
        [int(value) for value in row['expected_cp'].split(';') if value],
        dtype=np.int64,
    )
    np.testing.assert_array_equal(result.cp_set, expected)
    assert result.data.shape == data.shape
    assert not result.data.flags.writeable
    # Rank is transformed before entering the native mean family, matching
    # the public R result metadata.  Other fixture names are native families.
    expected_family = 'mean' if row['family'] == 'rank' else row['family']
    assert result.family == expected_family


@pytest.mark.parametrize(
    'case_id',
    tuple(_SHARED_GENERIC_DETECTOR_ROWS),
    ids=tuple(_SHARED_GENERIC_DETECTOR_ROWS),
)
def test_shared_fixture_wrapper_matches_generic_dispatch(case_id):
    """Generic R families and their wrappers share the fixture contract."""
    row = _SHARED_GENERIC_DETECTOR_ROWS[case_id]
    data = _shared_fixture_data(row['data_file'])
    wrapped = _run_shared_fixture(row, data)
    order = _shared_fixture_order(row['order'])
    beta = (
        row['beta']
        if row['beta'] in {'BIC', 'MBIC', 'MDL'}
        else float(row['beta'])
    )
    generic_kwargs = {
        'data': data,
        'family': row['family'],
        'beta': beta,
        'cost_adjustment': row['cost_adjustment'],
        'trim': float(row['trim']),
        'vanilla_percentage': float(row['vanilla_percentage']),
    }
    if row['family'] in {'arima', 'var'}:
        generic_kwargs['order'] = order
    elif row['family'] in {
        'lm', 'mean', 'variance', 'meanvariance', 'exponential'
    }:
        generic_kwargs['order'] = order
    generic = detect(**generic_kwargs)
    np.testing.assert_array_equal(wrapped.cp_set, generic.cp_set)


@pytest.mark.parametrize(
    'case_id',
    tuple(_SHARED_VARIANCE_ROWS),
    ids=tuple(_SHARED_VARIANCE_ROWS),
)
def test_shared_fixture_variance_estimators_are_deterministic(case_id):
    """Variance helpers consume manifest-selected data and expectations."""
    row = _SHARED_VARIANCE_ROWS[case_id]
    data = _shared_fixture_data(row['data_file'])
    expected = float(row['expected_value'])
    tolerance = float(row['tolerance'])

    if row['operation'] == 'estimate_variance_arma':
        order = _shared_fixture_order(row['order'])
        data = data[:, 0]
        direct = fastcpd_pkg.estimate_variance_arma(
            data, p=order[0], q=order[1]
        )
        generic = fastcpd_pkg.estimate_variance(
            data, family=row['family'], p=order[0], q=order[1]
        )
        max_order = order[0] * order[1]
        assert len(direct.table) == max_order
        assert [item['model'] for item in direct.table] == [
            f'AR({index})' for index in range(1, max_order + 1)
        ]
        direct_value = direct.sigma2_bic
        generic_value = generic.sigma2_bic
    else:
        family = row['family']
        if family == 'mean':
            direct = fastcpd_pkg.estimate_variance_mean(data)
        elif family == 'median':
            data = data[:, 0]
            direct = fastcpd_pkg.estimate_variance_median(data)
        elif family == 'lm':
            direct = fastcpd_pkg.estimate_variance_linear_regression(data)
        else:
            raise AssertionError(
                f'No variance fixture registered for {family!r}'
            )
        generic = fastcpd_pkg.estimate_variance(data, family=family)
        direct_value = direct
        generic_value = generic

    np.testing.assert_allclose(
        direct_value, expected, rtol=tolerance, atol=tolerance
    )
    np.testing.assert_allclose(
        generic_value, direct_value, rtol=tolerance, atol=tolerance
    )


def test_shared_fixture_files_have_no_missing_or_extra_manifest_rows():
    """Catch accidental fixture renames before the language suites diverge."""
    listed = {row['data_file'] for row in _SHARED_ROWS.values()}
    csv_files = {path.name for path in _SHARED_FIXTURE_DIR.glob('*.csv')}
    assert listed == csv_files


def test_shared_fixture_generator_matches_committed_content():
    """The generator reproduces fixture content across Git checkout EOLs."""
    completed = subprocess.run(
        [
            sys.executable,
            str(_SHARED_FIXTURE_DIR / 'generate_fixtures.py'),
            '--check',
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


if __name__ == "__main__":
    unittest.main()
