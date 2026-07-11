"""
Perform change point detection using fastcpd.
"""

import collections
import numpy
import fastcpd.variance_estimation
from fastcpd.interface import fastcpd_impl

# Families dispatched to the C++ Python binding (NO_RCPP mode).
#
# PELT families: mean, variance, meanvariance, exponential, mgaussian, garch.
# SEN families: lasso, gaussian/lm, binomial, poisson, arma, ma.
# 'arima' is accepted by detect() and pre-differenced before C++ dispatch.
_SUPPORTED_FAMILIES = frozenset({
    'mean', 'variance', 'meanvariance', 'exponential', 'mgaussian', 'lasso',
    'garch', 'gaussian', 'binomial', 'poisson', 'arma', 'ma', 'arima',
})

# Map R-style alias names to the internal C++ family string.
_FAMILY_ALIASES = {
    'var': 'mgaussian',
    'lm':  'gaussian',
}

class CpdResult(collections.namedtuple(
    'CpdResult',
    ['cp_set', 'raw_cp_set', 'cost_values', 'residuals', 'thetas'],
)):
    """Result object returned by detect() when cp_only=False.

    Fields:
        cp_set: Change-point indices (1-based, matching the R package).
        raw_cp_set: Raw change-point indices before boundary trimming.
        cost_values: Segment cost values.
        residuals: Nested list of shape (n_obs, n_response).
        thetas: Nested list of shape (n_params, n_segments); column j holds
            the estimated parameters for segment j.
    """

    __slots__ = ()

    def confint(self, *args, **kwargs):
        """Construct confidence intervals for this result.

        The Python result keeps fit output compact, so pass the original
        ``data`` and ``family`` arguments when calling this method.
        """
        from fastcpd.confidence import confint
        return confint(self, *args, **kwargs)


def mean(data, **kwargs):
    """Find change points efficiently in mean change models.

    Args:
        data: Univariate or multivariate data for mean change detection.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='mean', **kwargs)


def exponential(data, **kwargs):
    """Find change points efficiently in exponentially distributed data.

    Args:
        data: Univariate data where each observation is exponentially
            distributed; the rate parameter is allowed to change.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='exponential', **kwargs)


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


def lasso(data, **kwargs):
    """Find change points efficiently in LASSO regression models.

    Args:
        data: Data where the first column is the response.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='lasso', **kwargs)


def garch(data, order=(1, 1), **kwargs):
    """Find change points in GARCH(p, q) models.

    Args:
        data: Univariate time series, shape (n,) or (n, 1).
        order: Tuple (p, q) — GARCH and ARCH orders.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='garch', order=order, **kwargs)


def lm(data, **kwargs):
    """Find change points in ordinary linear regression models.

    Args:
        data: Array where column 0 is the response and the remaining columns
            are predictors, shape (n, p+1).
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='gaussian', **kwargs)


def binomial(data, **kwargs):
    """Find change points in logistic regression models.

    Args:
        data: Array where column 0 is the binary response and the remaining
            columns are predictors, shape (n, p+1).
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='binomial', **kwargs)


def poisson(data, **kwargs):
    """Find change points in Poisson regression models.

    Args:
        data: Array where column 0 is the count response and the remaining
            columns are predictors, shape (n, p+1).
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='poisson', **kwargs)


def arma(data, order=(1, 0), **kwargs):
    """Find change points in ARMA(p, q) models.

    When order[0] == 0 (pure MA), routes to the MA family automatically.

    Args:
        data: Univariate time series, shape (n,) or (n, 1).
        order: Tuple (p, q) — AR and MA orders.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='arma', order=order, **kwargs)


def ar(data, order=1, **kwargs):
    """Find change points in AR(p) models (pure autoregressive).

    Args:
        data: Univariate time series, shape (n,) or (n, 1).
        order: AR order p.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='arma', order=(order, 0), **kwargs)


def arima(data, order=(1, 1, 0), **kwargs):
    """Find change points in ARIMA(p, d, q) models.

    The integration order ``d`` is handled in Python by pre-differencing the
    series ``d`` times before running ARMA(p, q) change-point detection on the
    differenced series. This is equivalent to R's ``fastcpd_arima()`` for the
    common case ``d ≤ 2`` and matches its change-point indices exactly for
    ``d = 1`` (the most common case).

    Args:
        data: Univariate time series, shape (n,) or (n, 1).
        order: Tuple (p, d, q) — AR order, integration order, MA order.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.

    Note:
        For ``d ≥ 2`` the returned change-point indices correspond to the
        differenced-series positions. For ``d = 0`` this is identical to
        ``arma(data, order=(p, q))``.
    """
    return detect(data=data, family='arima', order=order, **kwargs)


def confint(result, **kwargs):
    """Construct confidence intervals for a ``CpdResult``.

    This is the Python analogue of R's ``confint(result, ...)`` API. The
    result object also exposes ``result.confint(...)``.
    """
    from fastcpd.confidence import confint as _confint
    return _confint(result, **kwargs)


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
    show_progress: bool = False,
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
        family: One of 'mean', 'variance', 'meanvariance', 'exponential',
            'mgaussian' / 'var' (alias), 'lasso', 'garch', 'gaussian' /
            'lm' (alias), 'binomial', 'poisson', 'arma', 'ma',
            'arima' (pre-differenced then routed to arma/ma).
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
        vanilla_percentage: Fraction of observations evaluated with pure PELT
            (no gradient update). 1.0 runs full PELT; 0.0 runs full SEN.
        warm_start: If True, use previous segment parameters as initial
            values.
        show_progress: If True, display a tqdm-format progress bar on stderr
            showing PELT timestep progress. Same format as Python tqdm default:
            ``42%|████████████          | 42/100 [00:05<00:07, 8.33it/s]``.
            Implemented in C++; no tqdm package required.

    Returns:
        When cp_only=True: a list of change-point indices (1-based).
        When cp_only=False: a CpdResult namedtuple.
    """
    if data is None:
        raise ValueError("data must be provided")
    data = numpy.asarray(data, dtype=float)
    if data.ndim == 1:
        data = data.reshape(-1, 1)

    if cost is not None or cost_gradient is not None or cost_hessian is not None:
        raise NotImplementedError(
            "Custom cost functions (cost, cost_gradient, cost_hessian) are not "
            "supported by the Python binding."
        )

    family = family.lower() if family is not None else 'custom'
    # Apply aliases (e.g. 'var' → 'mgaussian', 'lm' → 'gaussian').
    family = _FAMILY_ALIASES.get(family, family)

    # ARIMA(p, d, q): pre-difference the series d times, then route to arma/ma.
    # numpy.diff does not change the number of rows (it shortens by d rows),
    # so returned change-point indices are in the original-series index space.
    if family == 'arima':
        arima_order = list(order) if hasattr(order, '__len__') else [int(order), 0, 0]
        while len(arima_order) < 3:
            arima_order.append(0)
        p_ar, d_int, q_ma = int(arima_order[0]), int(arima_order[1]), int(arima_order[2])
        if d_int < 0:
            raise ValueError(f"ARIMA integration order d must be >= 0, got {d_int}")
        if d_int > 0:
            data = numpy.diff(data, n=d_int, axis=0)
        order = (p_ar, q_ma)
        family = 'arma' if p_ar > 0 else 'ma'

    # Route pure-MA models: arma with p=0 uses the MA family.
    if family == 'arma' and hasattr(order, '__len__') and order[0] == 0:
        family = 'ma'

    # Pure AR(p): route through lag-structured gaussian instead of ARMA.
    # The ARMA C++ path requires q > 0 in the NO_RCPP (Python) build;
    # for q == 0, OLS on lagged data gives the exact conditional MLE.
    _ar_offset = 0
    if (family == 'arma' and hasattr(order, '__len__') and
            len(order) >= 2 and int(order[0]) > 0 and int(order[1]) == 0):
        p_ar = int(order[0])
        n_rows = data.shape[0]
        if p_ar < n_rows:
            y = data[p_ar:, 0:1]
            lags = numpy.column_stack(
                [data[p_ar - j - 1:n_rows - j - 1, 0] for j in range(p_ar)])
            data = numpy.column_stack([y, lags])
            family = 'gaussian'
            _ar_offset = p_ar

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
        bool(show_progress),
    )

    if cp_only:
        return [cp + _ar_offset for cp in result['cp_set']]

    return CpdResult(
        cp_set=[cp + _ar_offset for cp in result['cp_set']],
        raw_cp_set=[cp + _ar_offset for cp in result['raw_cp_set']],
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
