import unittest

import numpy as np
from fastcpd.segmentation import exponential, lasso, mean, meanvariance, var, variance
from numpy import concatenate
from numpy.random import exponential as rexp
from numpy.random import multivariate_normal, randn, seed


class TestBasic(unittest.TestCase):

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


if __name__ == "__main__":
    unittest.main()
