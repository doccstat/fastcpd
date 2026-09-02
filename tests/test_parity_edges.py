"""Focused tests for cross-language covariance and traceback contracts."""

import numpy as np
import pytest

from fastcpd import detect_lm, detect_mean
from fastcpd.interface import fastcpd_impl as native_fastcpd_impl
from fastcpd.variance_estimation import estimate_variance_linear_regression


def test_raw_cp_set_keeps_boundary_removed_by_trim():
    """The native traceback remains available when trim removes a CP."""
    data = np.r_[np.zeros(1), np.full(9, 5.0)]
    result = detect_mean(
        data,
        beta=0.0,
        cost_adjustment="BIC",
        trim=0.1,
        variance_estimation=np.eye(1),
        cp_only=True,
    )

    np.testing.assert_array_equal(result.cp_set, np.empty(0, dtype=np.int64))
    np.testing.assert_array_equal(result.raw_cp_set, np.array([1], dtype=np.int64))


def test_raw_cp_traceback_is_ascending_and_matches_full_result():
    """The O(k) native traceback preserves every result-facing boundary."""
    data = np.r_[np.zeros(10), np.full(10, 5.0), np.full(10, -3.0)]
    options = {
        "beta": 5.0,
        "cost_adjustment": "BIC",
        "trim": 0.0,
        "variance_estimation": np.eye(1),
    }
    full = detect_mean(data, cp_only=False, **options)
    cp_only = detect_mean(data, cp_only=True, **options)

    expected = np.array([10, 20], dtype=np.int64)
    np.testing.assert_array_equal(full.cp_set, expected)
    np.testing.assert_array_equal(full.raw_cp_set, expected)
    np.testing.assert_array_equal(cp_only.cp_set, full.cp_set)
    np.testing.assert_array_equal(cp_only.raw_cp_set, full.raw_cp_set)
    assert np.all(np.diff(full.raw_cp_set) > 0)
    assert full.cost_values.shape == (3,)
    assert full.residuals.shape == (30, 1)
    assert full.thetas.shape == (1, 3)


def test_variance_lm_all_failed_blocks_returns_nan():
    """Singular regression blocks are omitted, as in R variance.lm()."""
    data = np.column_stack([
        np.arange(1.0, 6.0),
        np.ones(5),
        np.ones(5),
    ])
    estimate = estimate_variance_linear_regression(data)
    assert np.isnan(estimate)


def test_variance_lm_zero_iqr_matches_r_nan_threshold():
    """R's default Inf-IQR filter returns NaN for identical zero estimates."""
    predictor = np.arange(1.0, 12.0)
    data = np.column_stack([2.0 * predictor, predictor])
    assert np.isnan(estimate_variance_linear_regression(data))


def test_singular_explicit_gaussian_covariance_uses_rcond_fallback():
    """A singular scalar covariance produces the 1e-10 fallback result."""
    x = np.linspace(-1.0, 1.0, 40)
    y = np.sin(x)
    data = np.column_stack([y, x])
    options = {
        "beta": 1.0,
        "cost_adjustment": "BIC",
        "vanilla_percentage": 1.0,
        "cp_only": False,
    }
    singular = detect_lm(
        data, variance_estimation=np.array([[0.0]]), **options
    )
    fallback = detect_lm(
        data, variance_estimation=np.array([[1e-10]]), **options
    )
    for field in (
        "cp_set", "raw_cp_set", "cost_values", "residuals", "thetas"
    ):
        np.testing.assert_allclose(
            getattr(singular, field), getattr(fallback, field), equal_nan=True
        )


def test_character_beta_rejects_multielement_gaussian_covariance_like_r():
    """R's scalar beta rescaling rejects a 2x2 covariance for lm fits."""
    data = np.column_stack([
        np.linspace(-1.0, 1.0, 40),
        np.sin(np.linspace(-1.0, 1.0, 40)),
    ])
    with np.testing.assert_raises_regex(ValueError, 'must be scalar'):
        detect_lm(
            data,
            beta='BIC',
            variance_estimation=np.eye(2),
            cp_only=True,
        )


def test_native_character_beta_rejects_multielement_gaussian_covariance():
    """The standalone core enforces the scalar Gaussian penalty contract."""
    x = np.linspace(-1.0, 1.0, 40)
    data = np.column_stack([np.sin(x), x])
    with pytest.raises(ValueError, match="must be scalar"):
        native_fastcpd_impl(
            "BIC", "BIC", True, data, 1e-10, "gaussian",
            np.array([1.0]), np.empty(0), 0.0, np.zeros(3), 1, 1,
            float("nan"), 2, 0.0, np.empty(0), 1.0, np.eye(2), False,
            False,
        )


def test_nonfinite_covariance_is_rejected_instead_of_replaced():
    """R errors on nonfinite covariance input rather than changing the fit."""
    x = np.linspace(-1.0, 1.0, 40)
    with pytest.raises(ValueError, match="finite"):
        detect_lm(
            np.column_stack([np.sin(x), x]),
            beta=1.0,
            variance_estimation=np.array([[np.inf]]),
            cp_only=True,
        )
