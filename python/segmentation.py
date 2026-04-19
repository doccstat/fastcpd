"""
Perform change point detection using fastcpd.
"""

import collections
import numpy
import fastcpd.variance_estimation
from fastcpd.interface import fastcpd_impl

# Families supported by the C++ Python binding (NO_RCPP mode).
#
# PELT families (fully supported): mean, variance, meanvariance, mgaussian.
# SEN families (OLS-initialised, no glmnet/arima): lasso.
#
# Not yet supported (require R packages at runtime):
#   arma, arima, ma  — need stats::arima()
#   garch            — needs tseries::garch() + Fortran
#   gaussian/lm, binomial, poisson — need glmnet / fastglm (Rcpp + Eigen)
_SUPPORTED_FAMILIES = frozenset({
    'mean', 'variance', 'meanvariance', 'mgaussian', 'lasso',
})

# Map R-style alias names to the internal C++ family string.
_FAMILY_ALIASES = {
    'var': 'mgaussian',
}

# Result object returned by detect() when cp_only=False.
# Fields:
#   cp_set      – list of change-point indices (1-based, matching R package)
#   raw_cp_set  – raw change-point indices before boundary trimming
#   cost_values – list of segment cost values (one per segment)
#   residuals   – nested list of shape (n_obs, n_response)
#   thetas      – nested list of shape (n_params, n_segments); column j holds
#                 the estimated parameters for segment j
CpdResult = collections.namedtuple(
    'CpdResult',
    ['cp_set', 'raw_cp_set', 'cost_values', 'residuals', 'thetas'],
)


def mean(data, **kwargs):
    """Find change points efficiently in mean change models.

    Args:
        data: Univariate or multivariate data for mean change detection.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='mean', **kwargs)


def variance(data, **kwargs):
    """Find change points efficiently in variance change models.

    Args:
        data: Univariate or multivariate data for variance change detection.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='variance', **kwargs)


def meanvariance(data, **kwargs):
    """Find change points efficiently in mean and/or variance change models.

    Args:
        data: Univariate or multivariate data for mean and/or variance change
            detection.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='meanvariance', **kwargs)


def var(data, order, **kwargs):
    """Find change points efficiently in VAR (vector autoregression) models.

    Args:
        data: Multivariate time series data, shape (n, q+p) where q columns
            are responses and p columns are lagged predictors.
        order: Number of lagged predictors per response (p).
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='mgaussian', **kwargs)


def arma(data, order=(0, 0), **kwargs):
    """Find change points efficiently in ARMA models.

    Args:
        data: Univariate time series data.
        order: Tuple (p, q) for AR and MA orders.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='arma', order=order, **kwargs)


def lasso(data, **kwargs):
    """Find change points efficiently in LASSO regression models.

    Args:
        data: Data where the first column is the response.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='lasso', **kwargs)


def ma(data, order=0, **kwargs):
    """Find change points efficiently in MA (moving average) models.

    Args:
        data: Univariate time series data.
        order: MA order q.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='ma', order=(0, order), **kwargs)


def detect(
    formula: str = 'y ~ . - 1',
    data: numpy.ndarray = None,
    beta='MBIC',
    cost_adjustment: str = 'MBIC',
    family: str = None,
    cost=None,
    cost_gradient=None,
    cost_hessian=None,
    line_search=(1,),
    lower=None,
    upper=None,
    pruning_coef=None,
    segment_count: int = 10,
    trim: float = 0.05,
    momentum_coef: float = 0.0,
    multiple_epochs=lambda x: 0,
    epsilon: float = 1e-10,
    order=(0, 0, 0),
    p: int = None,
    p_response: int = 0,
    variance_estimation=None,
    cp_only: bool = False,
    vanilla_percentage: float = 1.0,
    warm_start: bool = False,
    **kwargs
):
    r"""Find change points efficiently.

    Args:
        formula: A formula string (unused; present for API parity with R).
        data: A NumPy array of shape (n, d) containing the data.
        beta: Penalty criterion. One of 'BIC', 'MBIC', 'MDL', or a float.
            The numeric value of the penalty is computed by C++ when a string
            is supplied, using the same formulae as the R package.
        cost_adjustment: One of 'BIC', 'MBIC', 'MDL'.
        family: One of 'mean', 'variance', 'meanvariance', 'mgaussian',
            'var' (alias for 'mgaussian'), 'lasso'.
        line_search: Values for line search step sizes.
        lower: Lower bound for parameters after each update.
        upper: Upper bound for parameters after each update.
        pruning_coef: Base pruning coefficient for the PELT algorithm.
            ``None`` (default) lets C++ compute the appropriate value
            automatically based on ``cost_adjustment`` and ``family``,
            matching R's ``get_pruning_coef()`` behaviour.
        segment_count: Initial guess for number of segments.
        trim: Trimming proportion for boundary change points.
        momentum_coef: Momentum coefficient for parameter updates.
        epsilon: Epsilon for numerical stability.
        order: Order for ARMA/MA models as tuple (ar_order, ma_order).
        p: Number of model parameters.  ``None`` (or 0) triggers automatic
            inference from ``family`` and the data dimensions in C++,
            matching the R package's per-family formulas.
        p_response: Number of response columns (mgaussian only).
        variance_estimation: Pre-specified variance/covariance matrix.
            When not supplied, estimated automatically (Rice estimator for
            mean/mgaussian, identity otherwise).
        cp_only: If True, return only change-point indices (list of floats).
            If False, return a CpdResult namedtuple with cp_set, raw_cp_set,
            cost_values, residuals, and thetas.
        vanilla_percentage: Fraction of data to run in pure PELT mode.
            Currently fixed at 1.0 (full PELT) for the Python binding.
        warm_start: If True, use previous segment parameters as initial
            values.

    Returns:
        When cp_only=True: a list of change-point indices (1-based).
        When cp_only=False: a CpdResult namedtuple.
    """
    if data is None:
        raise ValueError("data must be provided")
    data = numpy.asarray(data, dtype=float)
    if data.ndim == 1:
        data = data.reshape(-1, 1)

    family = family.lower() if family is not None else 'custom'
    # Apply aliases (e.g. 'var' → 'mgaussian').
    family = _FAMILY_ALIASES.get(family, family)

    if family not in _SUPPORTED_FAMILIES:
        raise ValueError(
            f"Family '{family}' is not supported by the Python binding. "
            f"Supported families: {sorted(_SUPPORTED_FAMILIES)}."
        )
    if cost_adjustment not in ('BIC', 'MBIC', 'MDL'):
        raise ValueError(
            f"cost_adjustment must be 'BIC', 'MBIC', or 'MDL', "
            f"got {cost_adjustment!r}"
        )

    # Variance estimation (NumPy; kept in Python because it uses array ops).
    ve = _estimate_variance(data, family, p_response, variance_estimation)

    # p=None or p=0 → C++ auto-infers from family + data shape.
    p_int = int(p) if (p is not None and p > 0) else 0

    # pruning_coef=None → C++ auto-computes (NaN is the sentinel).
    pruning_float = (
        float('nan') if pruning_coef is None else float(pruning_coef)
    )

    result = fastcpd_impl(
        beta,                   # str or float – C++ handles both
        cost_adjustment,
        bool(cp_only),
        data.tolist(),
        float(epsilon),
        family,
        list(line_search),
        list(lower) if lower is not None else [],
        float(momentum_coef),
        list(order) if hasattr(order, '__len__') else [float(order)],
        p_int,
        int(p_response),
        pruning_float,          # NaN → auto-compute in C++
        int(segment_count),
        float(trim),
        list(upper) if upper is not None else [],
        float(vanilla_percentage),
        ve.tolist(),
        bool(warm_start),
    )

    if cp_only:
        return result['cp_set']

    return CpdResult(
        cp_set=result['cp_set'],
        raw_cp_set=result['raw_cp_set'],
        cost_values=result['cost_values'],
        residuals=result['residuals'],
        thetas=result['thetas'],
    )


def _estimate_variance(data, family, p_response, variance_estimation):
    """Estimate the variance/covariance matrix for the given family."""
    if variance_estimation is not None:
        return numpy.asarray(variance_estimation, dtype=float)
    if family == 'mean':
        return fastcpd.variance_estimation.mean(data)
    if family == 'mgaussian':
        # Estimate Σ using Rice estimator on OLS residuals.
        # Using raw Y differences inflates Σ when predictors have high
        # variance, so we first partial out the predictor effect.
        d = data.shape[1]
        q = p_response if p_response > 0 else d
        y_cols = data[:, :q]
        x_cols = data[:, q:]
        if x_cols.shape[1] > 0:
            b_hat, _, _, _ = numpy.linalg.lstsq(x_cols, y_cols, rcond=None)
            resid = y_cols - x_cols @ b_hat
        else:
            resid = y_cols - y_cols.mean(axis=0)
        diffs = resid[1:] - resid[:-1]
        return numpy.mean(diffs[:, :, None] * diffs[:, None, :], axis=0) / 2
    # All other families: variance_estimate is not used by the C++ cost
    # function, so a 1×1 identity placeholder is sufficient.
    return numpy.eye(1)
