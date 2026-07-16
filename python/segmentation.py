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
# SEN families: lasso, gaussian/lm, binomial, poisson, quantile, arma, ma.
# 'arima' uses segment-local differencing in the shared C++ family.
# 'rank' and 'kernel'/'kcp' are Python-layer transforms routed to mean.
_SUPPORTED_FAMILIES = frozenset({
    'mean', 'variance', 'meanvariance', 'exponential', 'mgaussian', 'lasso',
    'garch', 'gaussian', 'binomial', 'poisson', 'quantile', 'arma', 'ma',
    'ar', 'arima', 'rank', 'kernel', 'kcp',
})

# Map R-style synonym names to the internal C++ family string.
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


def detect_mean(data, **kwargs):
    """Find change points efficiently in mean change models.

    Args:
        data: Univariate or multivariate data for mean change detection.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='mean', **kwargs)


def detect_exponential(data, **kwargs):
    """Find change points efficiently in exponentially distributed data.

    Args:
        data: Univariate data where each observation is exponentially
            distributed; the rate parameter is allowed to change.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='exponential', **kwargs)


def detect_variance(data, **kwargs):
    """Find change points efficiently in variance change models.

    Args:
        data: Univariate or multivariate data for variance change detection.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='variance', **kwargs)


def detect_meanvariance(data, **kwargs):
    """Find change points efficiently in mean and/or variance change models.

    Args:
        data: Univariate or multivariate data for mean and/or variance change
            detection.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='meanvariance', **kwargs)


def detect_var(data, order, **kwargs):
    """Find change points efficiently in VAR (vector autoregression) models.

    Args:
        data: Unlagged multivariate time series data, shape (n, q).
        order: Number of lagged predictors per response (p).
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='var', order=order, **kwargs)


def detect_lasso(data, **kwargs):
    """Find change points efficiently in LASSO regression models.

    Args:
        data: Data where the first column is the response.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='lasso', **kwargs)


def detect_garch(data, order=(1, 1), **kwargs):
    """Find change points in GARCH(p, q) models.

    Args:
        data: Univariate time series, shape (n,) or (n, 1).
        order: Tuple (p, q) — GARCH and ARCH orders.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='garch', order=order, **kwargs)


def detect_lm(data, **kwargs):
    """Find change points in ordinary linear regression models.

    Args:
        data: Array where column 0 is the response and the remaining columns
            are predictors, shape (n, p+1).
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='gaussian', **kwargs)


def detect_binomial(data, **kwargs):
    """Find change points in logistic regression models.

    Args:
        data: Array where column 0 is the binary response and the remaining
            columns are predictors, shape (n, p+1).
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='binomial', **kwargs)


def detect_poisson(data, **kwargs):
    """Find change points in Poisson regression models.

    Args:
        data: Array where column 0 is the count response and the remaining
            columns are predictors, shape (n, p+1).
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='poisson', **kwargs)


def detect_quantile(data, order=0.5, **kwargs):
    """Find change points in quantile regression models.

    Args:
        data: Array where column 0 is the response and the remaining columns
            are predictors, shape (n, p+1).
        order: Quantile level in (0, 1).
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    if not 0 < order < 1:
        raise ValueError(f"order must be in (0, 1), got {order!r}")
    return detect(data=data, family='quantile', order=(order,), **kwargs)


def detect_arma(data, order=(1, 0), **kwargs):
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


def detect_ar(data, order=1, **kwargs):
    """Find change points in AR(p) models (pure autoregressive).

    Args:
        data: Univariate time series, shape (n,) or (n, 1).
        order: AR order p.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    return detect(data=data, family='ar', order=order, **kwargs)


def detect_arima(data, order=(1, 1, 0), include_mean=False, **kwargs):
    """Find change points in ARIMA(p, d, q) models.

    The integration order ``d`` is applied independently inside every
    candidate segment by the shared native R/Python implementation. This
    avoids forming a difference across a proposed change-point boundary and
    keeps change-point indices in the original-series coordinate system.

    Args:
        data: Univariate time series, shape (n,) or (n, 1).
        order: Tuple (p, d, q) — AR order, integration order, MA order.
        include_mean: Must remain False. The unified likelihood is zero-mean.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.

    For ``d = 0`` this is identical to
    ``arma(data, order=(p, q))``.
    """
    return detect(
        data=data, family='arima', order=order,
        include_mean=include_mean, **kwargs
    )


def confint(result, **kwargs):
    """Construct confidence intervals for a ``CpdResult``.

    This is the Python analogue of R's ``confint(result, ...)`` API. The
    result object also exposes ``result.confint(...)``.
    """
    from fastcpd.confidence import confint as _confint
    return _confint(result, **kwargs)


def detect_rank(data, **kwargs):
    """Find change points using rank-transformed observations.

    Each column is replaced by its centered average rank, then mean-change
    detection is applied to the transformed data.
    """
    return detect(data=data, family='rank', **kwargs)


def detect_kernel(data, order=(100, 0), random_state=None, **kwargs):
    """Find distributional change points using random Fourier features.

    Args:
        data: Univariate or multivariate data.
        order: Tuple ``(n_features, bandwidth)``. A non-positive bandwidth
            uses the median heuristic on up to 1000 sampled observations.
        random_state: Optional NumPy random seed or generator.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A list of change-point indices, or a CpdResult when cp_only=False.
    """
    kwargs.setdefault('cost_adjustment', 'BIC')
    return detect(
        data=data, family='kernel', order=order, random_state=random_state,
        **kwargs
    )


def detect_kcp(data, order=(100, 0), random_state=None, **kwargs):
    """Find distributional change points using the KCP wrapper name."""
    return detect_kernel(
        data, order=order, random_state=random_state, **kwargs
    )


def detect_mean_variance(data, **kwargs):
    """Find change points in mean and/or variance change models."""
    return detect(data=data, family='meanvariance', **kwargs)


def detect_linear_regression(data, **kwargs):
    """Find change points in ordinary linear regression models."""
    return detect(data=data, family='gaussian', **kwargs)


def detect_logistic_regression(data, **kwargs):
    """Find change points in logistic regression models."""
    return detect(data=data, family='binomial', **kwargs)


def detect_poisson_regression(data, **kwargs):
    """Find change points in Poisson regression models."""
    return detect(data=data, family='poisson', **kwargs)


def detect_quantile_regression(data, order=0.5, **kwargs):
    """Find change points in quantile regression models."""
    return detect_quantile(data, order=order, **kwargs)


mean = detect_mean
exponential = detect_exponential
variance = detect_variance
meanvariance = detect_meanvariance
var = detect_var
lasso = detect_lasso
garch = detect_garch
lm = detect_lm
binomial = detect_binomial
poisson = detect_poisson
quantile = detect_quantile
arma = detect_arma
ar = detect_ar
arima = detect_arima
rank = detect_rank
kernel = detect_kernel
kcp = detect_kcp


def detect(
    formula: str = 'y ~ . - 1',
    data: numpy.ndarray = None,
    beta='MBIC',
    cost_adjustment: str = None,
    family: str = None,
    cost=None,
    cost_gradient=None,
    cost_hessian=None,
    line_search=(1,),
    lower=None,
    upper=None,
    pruning_coef=None,
    segment_count: int = 10,
    trim: float = 0.0,
    momentum_coef: float = 0.0,
    multiple_epochs=None,
    epsilon: float = 1e-10,
    order=(0, 0, 0),
    include_mean: bool = False,
    p: int = None,
    p_response: int = 0,
    variance_estimation=None,
    cp_only: bool = False,
    vanilla_percentage: float = 0.0,
    warm_start: bool = False,
    show_progress: bool = False,
    random_state=None,
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
            'mgaussian' / 'var' (synonym), 'lasso', 'garch', 'gaussian' /
            'lm' (synonym), 'binomial', 'poisson', 'arma', 'ma',
            'quantile', 'arima' (segment-local differencing),
            'rank', or 'kernel' / 'kcp' (random Fourier feature transform).
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
        multiple_epochs: Per-step epoch schedule. Custom schedules are not yet
            supported by the Python binding; leave this as ``None``.
        epsilon: Epsilon for numerical stability.
        order: Model order. ARMA uses ``(p, q)`` and ARIMA uses ``(p, d, q)``.
        include_mean: Must be False for ARIMA. The shared native likelihood
            is zero-mean in both R and Python.
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
    if multiple_epochs is not None:
        raise NotImplementedError(
            "Custom multiple_epochs schedules are not supported by the "
            "Python binding."
        )
    data = numpy.asarray(data, dtype=float)
    if data.ndim == 1:
        data = data.reshape(-1, 1)

    if cost is not None or cost_gradient is not None or cost_hessian is not None:
        raise NotImplementedError(
            "Custom cost functions (cost, cost_gradient, cost_hessian) are not "
            "supported by the Python binding."
        )

    cost_adjustment_missing = cost_adjustment is None
    if cost_adjustment is None:
        cost_adjustment = 'MBIC'

    family = family.lower() if family is not None else 'custom'
    index_offset = 0
    if family == 'var':
        var_order = _validate_var_order(order)
        legacy_columns = p_response * (var_order + 1)
        if p_response > 0 and data.shape[1] == legacy_columns:
            # Backward compatibility: older Python releases expected var()
            # input to contain [responses, lagged predictors].
            order = (var_order,)
            family = 'mgaussian'
        else:
            data, response_count = _var_regression_data(data, var_order)
            if p_response not in (0, response_count):
                raise ValueError(
                    "p_response must match the number of columns in VAR data"
                )
            p_response = response_count
            p = var_order * response_count ** 2
            order = (var_order,)
            index_offset = var_order
            family = 'mgaussian'
    elif family == 'ar':
        ar_order = _validate_ar_order(order)
        if data.shape[1] != 1:
            raise ValueError("AR data must be univariate")
        if data.shape[0] <= ar_order:
            raise ValueError(
                "AR order must be smaller than the number of rows"
            )
        n_rows = data.shape[0]
        response = data[ar_order:, 0:1]
        lags = numpy.column_stack([
            data[ar_order - lag:n_rows - lag, 0]
            for lag in range(1, ar_order + 1)
        ])
        data = numpy.column_stack([response, lags])
        family = 'gaussian'
        order = (ar_order,)
        index_offset = ar_order
    else:
        # Apply family-name synonyms (e.g. 'lm' → 'gaussian').
        family = _FAMILY_ALIASES.get(family, family)

    if family == 'rank':
        data = _rank_transform(data)
        family = 'mean'
        vanilla_percentage = 1.0

    if family in ('kernel', 'kcp'):
        original_dim = data.shape[1]
        data = _kernel_transform(data, order=order, random_state=random_state)
        family = 'mean'
        vanilla_percentage = 1.0
        if isinstance(beta, str):
            beta = (original_dim + 2) * numpy.log(data.shape[0]) / 2
        if variance_estimation is None:
            variance_estimation = numpy.eye(data.shape[1])
        if cost_adjustment_missing:
            cost_adjustment = 'BIC'

    # ARIMA(p, d, q): d=0 is exactly ARMA(p, q); d>0 stays as the
    # segment-local native ARIMA family so no cross-boundary difference is
    # introduced and returned indices remain in original-series coordinates.
    if family == 'arima':
        if include_mean:
            raise ValueError(
                "include_mean=True is not supported by the unified ARIMA "
                "likelihood; use the default False"
            )
        p_ar, d_int, q_ma = _validate_arima_order(order)
        if data.shape[1] != 1:
            raise ValueError("ARIMA data must be univariate")
        if data.shape[0] <= d_int:
            raise ValueError(
                "ARIMA integration order must be smaller than the number "
                "of rows"
            )
        if d_int == 0:
            order = (p_ar, q_ma)
            family = 'arma'
        else:
            order = (p_ar, d_int, q_ma)
    elif include_mean:
        raise ValueError("include_mean is only supported for ARIMA")

    if family == 'arma':
        if data.shape[1] != 1:
            raise ValueError("ARMA data must be univariate")
        order = _validate_arma_order(order)

    # Route pure-MA models: arma with p=0 uses the MA family.
    if family == 'arma' and hasattr(order, '__len__') and order[0] == 0:
        family = 'ma'

    # Pure AR(p) uses the same lag-structured Gaussian path as R's `ar` and
    # `arma(order = c(p, 0))` interfaces. The sequential ARMA update is not
    # defined for q=0.
    if (family == 'arma' and hasattr(order, '__len__') and
            len(order) >= 2 and int(order[0]) > 0 and int(order[1]) == 0):
        p_ar = int(order[0])
        if data.shape[1] != 1:
            raise ValueError("pure AR data must be univariate")
        n_rows = data.shape[0]
        if n_rows <= p_ar:
            raise ValueError(
                "AR order must be smaller than the number of rows"
            )
        response = data[p_ar:, 0:1]
        lags = numpy.column_stack([
            data[p_ar - lag:n_rows - lag, 0]
            for lag in range(1, p_ar + 1)
        ])
        data = numpy.column_stack([response, lags])
        family = 'gaussian'
        index_offset += p_ar

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
        return [cp + index_offset for cp in result['cp_set']]

    return CpdResult(
        cp_set=[cp + index_offset for cp in result['cp_set']],
        raw_cp_set=[cp + index_offset for cp in result['raw_cp_set']],
        cost_values=result['cost_values'],
        residuals=result['residuals'],
        thetas=result['thetas'],
    )


def _integer_order_values(order, expected_length, family):
    """Return a fixed-length tuple of non-negative integer order values."""
    if not hasattr(order, '__len__') or isinstance(order, str):
        raise ValueError(
            f"{family} order must contain {expected_length} integers"
        )
    values = list(order)
    if len(values) != expected_length:
        raise ValueError(
            f"{family} order must contain {expected_length} integers"
        )
    integers = []
    for value in values:
        value_float = float(value)
        value_int = int(value_float)
        if value_int < 0 or value_float != value_int:
            raise ValueError(
                f"{family} order must contain non-negative integers"
            )
        integers.append(value_int)
    return tuple(integers)


def _validate_ar_order(order):
    """Return a validated positive scalar AR order."""
    if hasattr(order, '__len__') and not isinstance(order, str):
        values = list(order)
        if len(values) != 1:
            raise ValueError("AR order must be a positive integer")
        order = values[0]
    value_float = float(order)
    value_int = int(value_float)
    if value_int <= 0 or value_float != value_int:
        raise ValueError("AR order must be a positive integer")
    return value_int


def _validate_arma_order(order):
    values = _integer_order_values(order, 2, "ARMA")
    if values == (0, 0):
        raise ValueError("ARMA order must contain at least one non-zero value")
    return values


def _validate_arima_order(order):
    values = _integer_order_values(order, 3, "ARIMA")
    if values == (0, 0, 0):
        raise ValueError("ARIMA order must contain at least one non-zero value")
    return values


def _validate_var_order(order):
    """Return a validated scalar VAR order."""
    if hasattr(order, '__len__') and not isinstance(order, str):
        order_values = list(order)
        if len(order_values) != 1:
            raise ValueError("VAR order must be a positive integer")
        order = order_values[0]
    order_float = float(order)
    order_int = int(order_float)
    if order_int <= 0 or order_float != order_int:
        raise ValueError("VAR order must be a positive integer")
    return order_int


def _var_regression_data(data, order):
    """Construct [responses, lagged predictors] for a VAR(p) series."""
    if data.ndim != 2 or data.shape[1] == 0:
        raise ValueError("VAR data must be a non-empty 2-D array")
    if data.shape[0] <= order:
        raise ValueError("VAR order must be smaller than the number of rows")

    responses = data[order:, :]
    predictors = numpy.column_stack([
        data[order - lag:data.shape[0] - lag, :]
        for lag in range(1, order + 1)
    ])
    return numpy.column_stack([responses, predictors]), data.shape[1]


def _rank_transform(data):
    data_matrix = numpy.asarray(data, dtype=float)
    if data_matrix.ndim == 1:
        data_matrix = data_matrix.reshape(-1, 1)
    ranks = numpy.column_stack([
        _average_ranks(data_matrix[:, column])
        for column in range(data_matrix.shape[1])
    ])
    return ranks - (data_matrix.shape[0] + 1) / 2


def _average_ranks(values):
    values = numpy.asarray(values, dtype=float)
    order = numpy.argsort(values, kind='mergesort')
    sorted_values = values[order]
    ranks = numpy.empty(values.shape[0], dtype=float)
    start = 0
    while start < values.shape[0]:
        end = start + 1
        while end < values.shape[0] and sorted_values[end] == sorted_values[start]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2
        start = end
    return ranks


def _kernel_transform(data, order=(100, 0), random_state=None):
    data_matrix = numpy.asarray(data, dtype=float)
    if data_matrix.ndim == 1:
        data_matrix = data_matrix.reshape(-1, 1)
    kernel_order = (
        list(order) if hasattr(order, '__len__') and not isinstance(order, str)
        else [order]
    )
    feature_count = (
        int(kernel_order[0]) if kernel_order and kernel_order[0] > 0 else 100
    )
    bandwidth = float(kernel_order[1]) if len(kernel_order) >= 2 else 0.0
    rng = _rng_from_random_state(random_state)
    if bandwidth <= 0:
        n_rows = data_matrix.shape[0]
        if n_rows > 1000:
            idx = rng.choice(n_rows, size=1000, replace=False)
            sampled = data_matrix[idx, :]
        else:
            sampled = data_matrix
        diffs = sampled[:, None, :] - sampled[None, :, :]
        squared_distances = numpy.sum(diffs * diffs, axis=2)
        positive_distances = squared_distances[squared_distances > 0]
        bandwidth = (
            numpy.sqrt(numpy.median(positive_distances) / 2)
            if positive_distances.size else 1.0
        )
    omega = rng.normal(
        loc=0.0, scale=1.0 / bandwidth,
        size=(data_matrix.shape[1], feature_count),
    )
    phase = rng.uniform(0.0, 2 * numpy.pi, size=feature_count)
    return numpy.sqrt(2.0 / feature_count) * numpy.cos(data_matrix @ omega + phase)


def _rng_from_random_state(random_state):
    if random_state is None:
        return numpy.random
    if isinstance(random_state, (numpy.random.Generator, numpy.random.RandomState)):
        return random_state
    return numpy.random.default_rng(random_state)


def _estimate_variance(data, family, p_response, variance_estimation):
    """Estimate the variance/covariance matrix for the given family."""
    if variance_estimation is not None:
        return numpy.asarray(variance_estimation, dtype=float)
    if family == 'mean':
        return fastcpd.variance_estimation.estimate_variance_mean(data)
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
