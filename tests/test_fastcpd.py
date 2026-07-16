import unittest

import fastcpd as fastcpd_pkg
import fastcpd.segmentation as segmentation
import numpy as np
from fastcpd.segmentation import (
    ar, arima, arma, binomial, detect_kernel, detect_mean, detect_quantile,
    detect_rank, exponential, garch, lasso, lm, mean, meanvariance, poisson,
    var, variance,
)
from numpy import concatenate
from numpy.random import exponential as rexp
from numpy.random import multivariate_normal, randn, seed


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
        self.assertEqual(result, [50.0])

    def test_unified_interface_aliases(self):
        self.assertIs(fastcpd_pkg.detect_mean, fastcpd_pkg.mean)
        self.assertIs(fastcpd_pkg.detect_kernel, fastcpd_pkg.kernel)
        self.assertIs(fastcpd_pkg.detect_kcp, fastcpd_pkg.kcp)
        self.assertIs(fastcpd_pkg.detect_lm, fastcpd_pkg.lm)
        self.assertIs(
            fastcpd_pkg.estimate_variance_mean,
            fastcpd_pkg.variance_estimation.mean,
        )
        for name in ('detect_time_series', 'detect_ts', 'time_series', 'ts'):
            self.assertFalse(hasattr(fastcpd_pkg, name))
            self.assertFalse(hasattr(segmentation, name))

        seed(17)
        data = concatenate((np.random.normal(0, 0.2, 40),
                            np.random.normal(3, 0.2, 40)))
        self.assertEqual(detect_mean(data).cp_set, mean(data).cp_set)
        self.assertEqual(
            detect_rank(data).cp_set,
            fastcpd_pkg.rank(data).cp_set,
        )
        self.assertEqual(
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
        self.assertEqual(
            fastcpd_pkg.detect_linear_regression(regression_data).cp_set,
            fastcpd_pkg.lm(regression_data).cp_set,
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

    def test_quantile_interface(self):
        seed(18)
        x = randn(120)
        y = concatenate((2 * x[:60], -2 * x[60:])) + 0.05 * randn(120)
        result = detect_quantile(np.column_stack([y, x]), order=0.5)
        self.assertGreater(len(result.cp_set), 0)
        self.assertAlmostEqual(result.cp_set[0], 60, delta=10)

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
            data=data, family='mean', method='profile', level=0.8, window=8)
        self.assertEqual(interval[0]['estimate'], result.cp_set[0])
        self.assertLessEqual(interval[0]['lower'], result.cp_set[0])
        self.assertGreaterEqual(interval[0]['upper'], result.cp_set[0])

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
        result = meanvariance(data)
        self.assertEqual(result.cp_set[0], 300)

    def test_var_mgaussian(self):
        # VAR(1) with 2 response cols: data = [y_t, y_{t-1}], shape (n-1, 4)
        seed(4)
        q = 2
        cov = [[1, 0], [0, 1]]
        y_raw = concatenate((
            multivariate_normal([0, 0], cov, 300),
            multivariate_normal([5, 5], cov, 300),
        ))
        data_mg = np.column_stack([y_raw[1:], y_raw[:-1]])
        result = var(data_mg, order=1, p_response=q)
        self.assertEqual(result.cp_set[0], 300)

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
        result = poisson(data)
        self.assertGreater(len(result.cp_set), 0)
        self.assertAlmostEqual(result.cp_set[0], 300, delta=20)

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
        result = arma(x, order=(1, 1))
        self.assertGreater(len(result.cp_set), 0)
        self.assertAlmostEqual(result.cp_set[0], 300, delta=30)

    def test_arima(self):
        # ARIMA(1,1,0): AR(1) process on first differences.
        # Build two I(1) segments with different AR(1) drift structures.
        seed(14)
        n = 600
        x = np.zeros(n)
        # Segment 1: I(1) with AR(1) φ=0.7 on differences
        for t in range(1, 300):
            x[t] = x[t - 1] + 0.7 * (x[t - 1] - x[t - 2] if t > 1 else 0) + randn()
        # Segment 2: I(1) with AR(1) φ=-0.7 on differences
        for t in range(300, n):
            x[t] = x[t - 1] - 0.7 * (x[t - 1] - x[t - 2] if t > 1 else 0) + randn()
        result = arima(x, order=(1, 1, 0))
        self.assertGreater(len(result.cp_set), 0)
        self.assertAlmostEqual(result.cp_set[0], 300, delta=30)

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
        self.assertEqual(r1.cp_set, r2.cp_set)


if __name__ == "__main__":
    unittest.main()
