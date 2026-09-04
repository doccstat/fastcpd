"""Focused regression tests for the cross-language Python API contract."""

import copy
import pickle
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from fastcpd import (
    confidence,
    detect_kernel,
    detect_mean,
    detect_meanvariance,
    detect_rank,
    lm,
    mv,
    var,
)
from fastcpd.segmentation import (
    CpdResult,
    _kernel_transform,
    _validate_kernel_order,
    arima,
)
from fastcpd import variance_estimation
from fastcpd.segmentation import detect


def test_variance_lm_validates_response_dimension_before_default_block_size():
    data = np.ones((12, 2))
    with pytest.raises(ValueError, match=r"d"):
        variance_estimation.estimate_variance_linear_regression(data, d=None)
    with pytest.raises(ValueError, match=r"d"):
        variance_estimation.estimate_variance_linear_regression(
            data, d=float("inf")
        )
    with pytest.raises(ValueError, match=r"d"):
        variance_estimation.estimate_variance_linear_regression(data, d=1.5)


def test_seeded_kernel_features_match_r_generation_order():
    fixture_dir = Path(__file__).parent / "fixtures" / "stochastic"
    data = np.loadtxt(
        fixture_dir / "kcp_seed_input.csv", delimiter=",", skiprows=1
    )
    expected = np.loadtxt(
        fixture_dir / "kcp_seed_features.csv", delimiter=",", skiprows=1
    )
    actual = _kernel_transform(data, order=(8, 1.25), random_state=7)
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=2e-15)


def test_identical_kernel_input_uses_shared_finite_bandwidth_fallback():
    result = detect(
        np.full(20, 3.0),
        family="kcp",
        order=(8, 0),
        random_state=7,
        beta=2,
        trim=0,
        cp_only=True,
    )
    np.testing.assert_array_equal(result.cp_set, np.empty(0, dtype=np.int64))


def test_seeded_mean_bootstrap_matches_r_interval():
    data = np.loadtxt(
        Path(__file__).parent / "fixtures" / "mean_step.csv",
        delimiter=",",
        skiprows=1,
    )
    result = detect_mean(
        data,
        beta=5,
        cost_adjustment="BIC",
        trim=0.1,
        vanilla_percentage=1,
    )
    interval = result.confint(
        parm="cp", method="bootstrap", B=12, level=0.8, random_state=17
    )
    assert interval == [{
        "parm": "cp",
        "index": 1,
        "estimate": 20,
        "lower": 20.0,
        "upper": 20.0,
        "detection_rate": 1.0,
        "level": 0.8,
        "method": "bootstrap",
        "bootstrap": "nonparametric",
    }]


def test_seeded_kcp_bootstrap_continues_the_r_stream():
    data = np.loadtxt(
        Path(__file__).parent / "fixtures" / "stochastic" /
        "kcp_seed_input.csv",
        delimiter=",",
        skiprows=1,
    )
    result = detect_kernel(
        data, order=(8, 1.25), random_state=7, beta=2, trim=0
    )
    interval = result.confint(
        parm="cp", method="bootstrap", B=12, level=0.8, random_state=19
    )
    assert interval == [{
        "parm": "cp",
        "index": 1,
        "estimate": 13,
        "lower": 11.0,
        "upper": 13.0,
        "detection_rate": 5 / 12,
        "level": 0.8,
        "method": "bootstrap",
        "bootstrap": "nonparametric",
    }]


def test_variance_lm_rejects_computationally_singular_blocks_like_r():
    # R's solve.default rejects these nearly collinear normal-equation
    # matrices because their reciprocal one-norm condition is below machine
    # epsilon.  NumPy's solve otherwise accepts them and can return a small
    # negative variance estimate.
    rng = np.random.default_rng(2)
    x = np.arange(1.0, 21.0)
    z = x + 1e-14 * rng.normal(size=x.size)
    data = np.column_stack([rng.normal(size=x.size), x, z])
    estimate = variance_estimation.estimate_variance_linear_regression(data)
    assert np.isnan(estimate)


def test_multivariate_variance_lm_all_failed_is_quiet_nan_matrix():
    """R quietly returns an all-NaN matrix when every block solve fails."""
    data = np.column_stack([
        np.arange(1.0, 8.0),
        np.arange(2.0, 16.0, 2.0),
        np.ones(7),
        np.ones(7),
    ])
    with np.errstate(all="raise"):
        estimate = variance_estimation.estimate_variance_linear_regression(
            data, d=2
        )
    assert estimate.shape == (2, 2)
    assert np.all(np.isnan(estimate))


def test_bootstrap_refit_removes_control_keyword_duplicates_and_supports_randomstate():
    calls = []

    class FakeSegmentation:
        @staticmethod
        def detect(*args, **kwargs):
            calls.append(("detect", kwargs))
            return SimpleNamespace(cp_set=np.empty(0, dtype=np.int64))

        @staticmethod
        def detect_kcp(*args, **kwargs):
            calls.append(("kcp", kwargs))
            return SimpleNamespace(cp_set=np.empty(0, dtype=np.int64))

    # Old/current fit metadata and explicit overrides may both contain these
    # controls.  They must be consumed exactly once by the bootstrap driver.
    confidence._bootstrap_refit(
        FakeSegmentation,
        np.zeros((8, 1)),
        "mean",
        {
            "family": "mean",
            "cp_only": False,
            "show_progress": True,
        },
        np.random.RandomState(7),
    )
    assert calls[-1][0] == "detect"
    assert calls[-1][1]["cp_only"] is True
    assert calls[-1][1]["show_progress"] is False

    confidence._bootstrap_refit(
        FakeSegmentation,
        np.zeros((8, 1)),
        "kcp",
        {
            "family": "kcp",
            "cp_only": False,
            "show_progress": True,
        },
        np.random.RandomState(7),
    )
    assert calls[-1][0] == "kcp"
    assert calls[-1][1]["cp_only"] is True
    assert "show_progress" not in calls[-1][1]
    assert isinstance(calls[-1][1]["random_state"], int)


def test_rank_bootstrap_refit_repeats_the_transform():
    """A rank wrapper result must not bootstrap raw values as mean data."""
    calls = []

    class FakeSegmentation:
        @staticmethod
        def detect_rank(*args, **kwargs):
            calls.append(("rank", args, kwargs))
            return SimpleNamespace(cp_set=np.empty(0, dtype=np.int64))

        @staticmethod
        def detect(*args, **kwargs):
            calls.append(("detect", args, kwargs))
            return SimpleNamespace(cp_set=np.empty(0, dtype=np.int64))

    confidence._bootstrap_refit(
        FakeSegmentation,
        np.zeros((8, 1)),
        "mean",
        {
            "family": "mean",
            "_wrapper_family": "rank",
            "_native_family": "mean",
        },
        np.random.default_rng(3),
    )
    assert calls[-1][0] == "rank"
    assert calls[-1][2]["cp_only"] is True


def test_rank_wald_uses_centered_ranks_like_r_result_data():
    raw = np.array([0.0, 100.0, 1.0, 101.0, 2.0, 102.0, 3.0, 103.0])
    result = detect_rank(
        raw,
        beta=1e6,
        cost_adjustment="BIC",
        trim=0,
        vanilla_percentage=1,
    )
    interval = result.confint(parm="theta", method="wald")
    centered_ranks = confidence._analysis_data(
        raw.reshape(-1, 1), result.fit_kwargs
    )[:, 0]
    expected_se = np.std(centered_ranks, ddof=1) / np.sqrt(raw.size)
    assert interval[0]["se"] == pytest.approx(expected_se)


def test_result_fit_options_are_immutable_pickleable_and_deepcopyable():
    result = detect_meanvariance(
        np.arange(12.0), beta=5.0, variance_estimation=np.eye(1)
    )
    with pytest.raises(TypeError):
        result.fit_kwargs["beta"] = 10.0

    restored = pickle.loads(pickle.dumps(result))
    duplicated = copy.deepcopy(result)
    for candidate in (restored, duplicated):
        assert set(candidate.fit_kwargs) == set(result.fit_kwargs)
        for key, expected in result.fit_kwargs.items():
            actual = candidate.fit_kwargs[key]
            if isinstance(expected, np.ndarray):
                np.testing.assert_array_equal(actual, expected)
                assert not actual.flags.writeable
            else:
                assert actual == expected
        np.testing.assert_array_equal(candidate.cp_set, result.cp_set)
        assert not candidate.data.flags.writeable


def test_manual_result_construction_still_copies_mutable_arrays():
    """The detector-only no-copy path must not weaken the public constructor."""
    source = np.arange(8.0)
    result = CpdResult(
        np.array([4]), np.array([4]), np.array([1.0]),
        np.zeros((8, 1)), np.array([[2.0]]), source,
    )
    source[0] = 999.0
    assert result.data[0, 0] == 0.0
    assert not result.data.flags.writeable


def test_legacy_five_field_result_can_supply_confint_context_explicitly():
    # The pre-metadata CpdResult had these five positional fields.  Keep that
    # construction valid while allowing confidence callers to provide the
    # original data/family/order explicitly.
    result = CpdResult(
        np.array([5]),
        np.array([5]),
        np.array([1.0]),
        np.zeros(10),
        np.array([[2.0]]),
    )
    assert result.data.shape == (0, 0)
    interval = confidence.confint(
        result,
        data=np.arange(10.0),
        family="mean",
        order=(0, 0, 0),
        method="profile",
        level=0.8,
        window=1,
    )
    assert interval[0]["estimate"] == 5


def test_cp_only_lagged_models_keep_original_data_dimensions():
    rng = np.random.default_rng(101)
    series = np.cumsum(rng.normal(size=40))
    arima_result = arima(series, order=(0, 1, 0), cp_only=True)
    assert arima_result.data.shape == (40, 1)
    assert arima_result.cost_values.shape == (0,)
    assert arima_result.residuals.shape == (0, 0)
    assert arima_result.thetas.shape == (0, 0)

    raw_var = rng.normal(size=(40, 2))
    var_result = var(raw_var, order=2, cp_only=True)
    assert var_result.data.shape == raw_var.shape
    assert var_result.cost_values.shape == (0,)
    assert var_result.residuals.shape == (0, 0)
    assert var_result.thetas.shape == (0, 0)


def test_multivariate_lm_wald_matches_flattened_theta_dimension():
    rng = np.random.default_rng(102)
    x = rng.normal(size=(80, 2))
    y = np.column_stack([
        x @ np.array([1.0, -0.5]),
        x @ np.array([0.25, 1.5]),
    ]) + 0.1 * rng.normal(size=(80, 2))
    result = lm(
        np.column_stack([y, x]),
        p_response=2,
        beta=1e6,
    )
    intervals = result.confint(parm="theta", method="wald")
    assert len(intervals) == result.thetas.size
    assert all(np.isfinite(row["se"]) for row in intervals)


def test_multivariate_residuals_keep_response_columns_and_r_flatten_order():
    """Python keeps the native n-by-q residual matrix.

    R's current wrapper applies ``matrix(result$residual)``, which flattens a
    multivariate matrix in column-major order.  The Python contract retains
    the more useful n-by-q shape; callers needing the R serialization can use
    ``result.residuals.ravel(order='F')`` and obtain exactly the same values.
    """
    rng = np.random.default_rng(103)
    x = rng.normal(size=(36, 2))
    y = np.column_stack([
        x @ np.array([1.0, -0.5]),
        x @ np.array([0.25, 1.5]),
    ]) + 0.05 * rng.normal(size=(36, 2))
    result = lm(np.column_stack([y, x]), p_response=2, beta=1e6)

    assert result.residuals.ndim == 2
    assert result.residuals.shape == (36, 2)
    np.testing.assert_array_equal(
        result.residuals.ravel(order='F'),
        np.concatenate((result.residuals[:, 0], result.residuals[:, 1])),
    )


def test_var_residuals_pad_each_response_in_original_coordinates():
    """Python keeps a complete n-by-q VAR residual matrix after lagging."""
    rng = np.random.default_rng(104)
    data = rng.normal(size=(36, 2))
    result = var(
        data,
        order=2,
        beta=1e6,
        variance_estimation=np.eye(2),
    )
    assert result.residuals.shape == data.shape
    assert np.all(np.isnan(result.residuals[:2, :]))
    assert np.all(np.isfinite(result.residuals[2:, :]))


def test_package_exports_mv_legacy_alias():
    assert mv is detect_meanvariance


def test_variance_lm_alias_is_direct_identity():
    """The R ``estimate_variance_lm`` spelling is a direct alias."""
    assert (
        variance_estimation.estimate_variance_lm
        is variance_estimation.estimate_variance_linear_regression
    )


def test_variance_aliases_are_direct_identities():
    """All underscore variance spellings remain portable aliases."""
    assert variance_estimation.variance_mean is (
        variance_estimation.estimate_variance_mean
    )
    assert variance_estimation.variance_median is (
        variance_estimation.estimate_variance_median
    )
    assert variance_estimation.variance_lm is (
        variance_estimation.estimate_variance_linear_regression
    )
    assert variance_estimation.variance_linear_regression is (
        variance_estimation.estimate_variance_linear_regression
    )
    assert variance_estimation.variance_arma is (
        variance_estimation.estimate_variance_arma
    )


def test_cost_adjustment_signature_default_is_public_mbic_string():
    """The internal omission sentinel compares like the R default."""
    import inspect

    default = inspect.signature(detect).parameters['cost_adjustment'].default
    assert default == 'MBIC'
    assert repr(default) == "'MBIC'"


def test_kcp_cost_adjustment_normalizes_omitted_and_none_to_bic():
    """KCP follows R's BIC adjustment for omitted and explicit NULL values."""
    data = np.column_stack([
        np.r_[np.zeros(12), np.ones(12)],
        np.r_[np.zeros(12), np.ones(12)],
    ])
    from fastcpd.segmentation import detect, detect_kcp

    omitted = detect_kcp(data, order=(5, 1), random_state=7, cp_only=True)
    explicit_none = detect_kcp(
        data, order=(5, 1), random_state=7, cost_adjustment=None,
        cp_only=True,
    )
    generic_none = detect(
        data, family='kcp', order=(5, 1), random_state=7,
        cost_adjustment=None, cp_only=True,
    )
    assert omitted.fit_kwargs['cost_adjustment'] == 'BIC'
    assert explicit_none.fit_kwargs['cost_adjustment'] == 'BIC'
    assert generic_none.fit_kwargs['cost_adjustment'] == 'BIC'


@pytest.mark.parametrize('family', ('rank', 'kernel'))
def test_generic_detect_rejects_wrapper_only_distribution_families(family):
    """R exposes rank/kernel through wrappers; generic detect accepts kcp."""
    with pytest.raises(ValueError, match="not supported"):
        detect(np.arange(12.0), family=family, cp_only=True)


def test_summary_reports_response_count_not_regression_column_count():
    predictor = np.arange(1.0, 31.0)
    scalar = lm(
        np.column_stack([2.0 * predictor + 0.1 * (predictor % 3), predictor]),
        beta=1e6,
    )
    assert scalar.summary()['n_response'] == 1

    predictors = np.column_stack([predictor, predictor ** 2 / 100])
    responses = np.column_stack([
        predictors @ np.array([1.0, -0.25]),
        predictors @ np.array([-0.5, 0.75]),
    ])
    multivariate = lm(
        np.column_stack([responses, predictors]),
        p_response=2,
        beta=1e6,
        variance_estimation=np.eye(2),
    )
    assert multivariate.summary()['n_response'] == 2


@pytest.mark.parametrize(
    ('order', 'expected'),
    [
        (None, (100, 0.0)),
        ((), (100, 0.0)),
        (3, (3, 0.0)),
        ((-1, 1), (100, 1.0)),
        ((0, 1), (100, 1.0)),
        ((1, -1), (1, -1.0)),
        ((3, 1, 999), (3, 1.0)),
    ],
)
def test_kcp_order_normalization_matches_r_defaults(order, expected):
    """Valid KCP order spellings retain R's defaults and fallback rules."""
    assert _validate_kernel_order(order) == expected


@pytest.mark.parametrize(
    'order',
    [
        ('bad', 1),
        (1, 'bad'),
        (0.5, 1),
        (float('nan'), 1),
        (float('inf'), 1),
        (3, float('nan')),
        (3, float('inf')),
    ],
)
def test_kcp_order_rejects_malformed_numeric_values(order):
    """Malformed KCP orders fail before native allocation with ValueError."""
    with pytest.raises(ValueError, match='KCP'):
        _validate_kernel_order(order)


def test_kcp_transform_uses_validated_feature_count_and_is_finite():
    """The validated order controls a finite, reproducible feature matrix."""
    data = np.arange(12.0).reshape(-1, 1)
    transformed = _kernel_transform(data, order=(3, 1, 99), random_state=7)
    assert transformed.shape == (12, 3)
    assert np.all(np.isfinite(transformed))


def test_confidence_quantiles_match_r_double_precision_values():
    """Profile/Wald diagnostics use the same qnorm/qchisq values as R."""
    assert confidence._normal_quantile(0.975) == pytest.approx(
        1.9599639845400534, abs=1e-15
    )
    assert confidence._chisq1(0.8) / 2 == pytest.approx(
        0.8211872075749078, abs=2e-15
    )


def test_generic_detect_accepts_positional_array_with_keyword_family():
    """The common Python ``detect(array, family=...)`` spelling is valid."""
    data = np.r_[np.zeros(12), np.ones(12)]
    result = detect_meanvariance(
        data, beta=1e6, variance_estimation=1, cp_only=True
    )
    # ``detect`` is imported lazily here to keep the module's public imports
    # focused on the aliases under test.
    from fastcpd.segmentation import detect
    positional = detect(
        data, family='meanvariance', beta=1e6,
        variance_estimation=1, cp_only=True,
    )
    np.testing.assert_array_equal(positional.cp_set, result.cp_set)
