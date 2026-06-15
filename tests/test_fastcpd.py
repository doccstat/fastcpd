import unittest

from fastcpd.segmentation import exponential, mean
from numpy import concatenate
from numpy.random import exponential as rexp
from numpy.random import multivariate_normal, seed


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


if __name__ == "__main__":
    unittest.main()
