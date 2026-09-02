"""
Perform change point detection using fastcpd.
"""

import dataclasses
from typing import Mapping

import numpy
from fastcpd import variance_estimation as _variance_estimation
from fastcpd.interface import fastcpd_impl

# Families dispatched to the C++ Python binding (NO_RCPP mode).
#
# PELT families: mean, variance, meanvariance, exponential, mgaussian, garch.
# SEN families: lasso, gaussian/lm, binomial, poisson, quantile, arma, ma.
# 'arima' uses segment-local differencing in the shared C++ family.
# 'kcp' is a Python-layer transform routed to mean. Rank and ``kernel`` are
# wrapper-only spellings, matching R's generic ``detect()`` family contract.
_SUPPORTED_FAMILIES = frozenset({
    'mean', 'variance', 'meanvariance', 'exponential', 'mgaussian', 'lasso',
    'garch', 'gaussian', 'binomial', 'poisson', 'quantile', 'arma', 'ma',
    'ar', 'arima', 'var', 'lm', 'kcp', 'custom',
})

# Map R-style synonym names to the internal C++ family string.
_FAMILY_ALIASES = {
    'var': 'mgaussian',
    'lm':  'gaussian',
}


class FrozenMapping(Mapping):
    """Small immutable, pickleable mapping used for stored fit options."""

    __slots__ = ('_data',)

    def __init__(self, value=None):
        self._data = {
            key: _freeze_value(item) for key, item in dict(value or {}).items()
        }

    def __getitem__(self, key):
        return self._data[key]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def __repr__(self):
        return f'FrozenMapping({self._data!r})'

    def __reduce__(self):
        return (type(self), (self._data,))

    def __deepcopy__(self, memo):
        return self


class _UnsetCostAdjustment(str):
    """Sentinel that prints like the public MBIC default.

    R distinguishes an omitted ``cost_adjustment`` argument (MBIC) from an
    explicitly supplied ``NULL`` (BIC).  Python's normal ``None`` default
    cannot make that distinction, so the public signature uses this sentinel.
    Its representation is the documented default so help()/inspect remain
    understandable without exposing the implementation detail.
    """

    def __new__(cls):
        return str.__new__(cls, 'MBIC')

    def __repr__(self):  # pragma: no cover - cosmetic introspection helper
        return "'MBIC'"


_COST_ADJUSTMENT_UNSET = _UnsetCostAdjustment()


def _freeze_value(value):
    """Return an immutable snapshot suitable for ``CpdResult.fit_kwargs``."""
    if isinstance(value, numpy.ndarray):
        copied = numpy.array(value, copy=True)
        copied.setflags(write=False)
        return copied
    if isinstance(value, Mapping):
        return FrozenMapping(value)
    if isinstance(value, list):
        return tuple(_freeze_value(item) for item in value)
    if isinstance(value, tuple):
        return tuple(_freeze_value(item) for item in value)
    if isinstance(value, set):
        return frozenset(_freeze_value(item) for item in value)
    return value


def _freeze_mapping(value):
    """Normalize an optional mapping to an immutable stored snapshot."""
    if value is None:
        return FrozenMapping()
    if not isinstance(value, Mapping):
        raise TypeError("fit_kwargs must be a mapping")
    return _freeze_value(value)


def _as_2d_data(data):
    """Convert array-like input to a C-contiguous two-dimensional matrix."""
    if data is None:
        raise ValueError("data must be provided")
    matrix = numpy.array(data, dtype=float, copy=True, order='C')
    if matrix.ndim == 1:
        matrix = matrix.reshape(-1, 1)
    if matrix.ndim != 2:
        raise ValueError("data must be one- or two-dimensional")
    if matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ValueError("data must be non-empty")
    if not numpy.all(numpy.isfinite(matrix)):
        raise ValueError("data must contain only finite numeric values")
    return matrix


def _as_univariate_series(data):
    """Flatten a time-series wrapper input using R's column-major order."""
    values = numpy.asarray(data, dtype=float)
    if values.ndim == 0 or values.ndim > 2:
        raise ValueError("time-series data must be one- or two-dimensional")
    series = values.ravel(order='F')
    if series.size == 0:
        raise ValueError("time-series data must be non-empty")
    if not numpy.all(numpy.isfinite(series)):
        raise ValueError("time-series data must contain only finite values")
    return numpy.ascontiguousarray(series)


def _with_metadata(result, data, *, family, order, fit_kwargs=None):
    """Reuse an immutable result while replacing wrapper-visible metadata."""
    return CpdResult(
        cp_set=result.cp_set,
        raw_cp_set=result.raw_cp_set,
        cost_values=result.cost_values,
        residuals=result.residuals,
        thetas=result.thetas,
        data=data,
        family=family,
        order=order,
        cp_only=result.cp_only,
        fit_kwargs=result.fit_kwargs if fit_kwargs is None else fit_kwargs,
        _copy_arrays=False,
    )


def _fit_kwargs_for_wrapper(family, order, kwargs):
    """Build a minimal immutable-friendly fit option map for confidence refits."""
    options = dict(kwargs)
    options['family'] = family
    options['order'] = tuple(order)
    return options


def _build_fit_kwargs(
    *, family, order, beta, cost_adjustment, line_search, lower, upper,
    pruning_coef, segment_count, trim, momentum_coef, epsilon, p,
    p_response, variance_estimation, vanilla_percentage, warm_start,
    show_progress, random_state=None,
):
    """Capture enough of a detection call for a faithful confidence refit."""
    options = {
        'family': family,
        'order': tuple(order),
        'beta': beta,
        'cost_adjustment': cost_adjustment,
        'line_search': line_search,
        'segment_count': segment_count,
        'trim': trim,
        'momentum_coef': momentum_coef,
        'epsilon': epsilon,
        'vanilla_percentage': vanilla_percentage,
        'warm_start': bool(warm_start),
        'show_progress': bool(show_progress),
    }
    if lower is not None:
        options['lower'] = lower
    if upper is not None:
        options['upper'] = upper
    if pruning_coef is not None:
        options['pruning_coef'] = pruning_coef
    if p is not None:
        options['p'] = p
    if p_response:
        options['p_response'] = p_response
    # Omitting an estimated covariance is important: bootstrap samples must
    # re-estimate it rather than reusing the original sample's value.
    if variance_estimation is not None:
        options['variance_estimation'] = variance_estimation
    if random_state is not None:
        # A NumPy Generator/RandomState is mutable and cannot be made truly
        # immutable by the recursive fit-options snapshot.  KCP bootstrap
        # refits draw a fresh scalar seed from their own RNG, so retaining the
        # live object here would only create an alias (and could mutate state
        # through ``result.fit_kwargs``).  Scalar seeds remain available for
        # introspection and legacy callers.
        if not isinstance(
            random_state, (numpy.random.Generator, numpy.random.RandomState)
        ):
            options['random_state'] = random_state
    return options


def _detect_kcp_wrapper(data, *, order, random_state, **kwargs):
    """Implement both Python KCP spellings on top of the R-compatible path."""
    data_matrix = _as_2d_data(data)
    public_order = _public_order(order)
    kwargs.setdefault('cost_adjustment', 'BIC')
    # R's KCP branch always uses the unadjusted BIC adjustment when callers
    # omit the argument or explicitly pass ``NULL``.  Normalize both Python
    # spellings before recording fit metadata so confidence refits and result
    # introspection describe the effective native option rather than leaking
    # the ``None`` sentinel through the wrapper boundary.
    if kwargs['cost_adjustment'] is None:
        kwargs['cost_adjustment'] = 'BIC'
    original_dim = data_matrix.shape[1]
    transformed = _kernel_transform(
        data_matrix, order=public_order, random_state=random_state
    )
    beta = kwargs.pop('beta', 'MBIC')
    beta_for_fit = beta
    if isinstance(beta, str):
        # R's KCP wrapper uses the original dimension in its MBIC penalty,
        # irrespective of the number of random features D.
        beta = (original_dim + 2) * numpy.log(data_matrix.shape[0]) / 2
    kwargs.setdefault('variance_estimation', numpy.eye(transformed.shape[1]))
    kwargs['vanilla_percentage'] = 1.0
    result = detect(
        data=transformed,
        family='mean',
        beta=beta,
        order=(0, 0, 0),
        **kwargs,
    )
    fit_kwargs = _build_fit_kwargs(
        family='kcp', order=public_order, beta=beta_for_fit,
        cost_adjustment=kwargs.get('cost_adjustment', 'BIC'),
        line_search=kwargs.get('line_search', (1,)),
        lower=kwargs.get('lower'), upper=kwargs.get('upper'),
        pruning_coef=kwargs.get('pruning_coef'),
        segment_count=kwargs.get('segment_count', 10),
        trim=kwargs.get('trim', 0.0),
        momentum_coef=kwargs.get('momentum_coef', 0.0),
        epsilon=kwargs.get('epsilon', 1e-10),
        p=kwargs.get('p'), p_response=kwargs.get('p_response', 0),
        variance_estimation=kwargs.get('variance_estimation'),
        vanilla_percentage=1.0, warm_start=kwargs.get('warm_start', False),
        show_progress=kwargs.get('show_progress', False),
        random_state=random_state,
    )
    return _with_metadata(
        result, data_matrix, family='kcp', order=public_order,
        fit_kwargs=fit_kwargs,
    )


@dataclasses.dataclass(frozen=True, slots=True, eq=False)
class CpdResult:
    """Stable result object returned by every ``detect()`` call.

    Fields:
        cp_set: Read-only int64 array of 1-based change-point indices.
        raw_cp_set: Read-only int64 array before boundary trimming.
        cost_values: Read-only float array of segment cost values.
        residuals: Read-only float array of shape (n_obs, n_response), aligned
            to the original input coordinates (leading model-lag rows are
            NaN, matching the R result).
        thetas: Read-only float array of shape (n_params, n_segments); column j
            holds
            the estimated parameters for segment j.
        data: Read-only copy of the original input data, always two-dimensional.
        family: Public family name requested by the caller.
        order: Public model order supplied by the caller.
        cp_only: Whether detailed native fit output was skipped.
        fit_kwargs: Read-only mapping of the original fit options. It is used
            by ``confint()`` to reproduce bootstrap refits faithfully.
    """

    cp_set: numpy.ndarray
    raw_cp_set: numpy.ndarray
    cost_values: numpy.ndarray
    residuals: numpy.ndarray
    thetas: numpy.ndarray
    # ``data``, ``family`` and ``order`` were added after the original
    # five-field result object.  Defaults preserve positional/keyword
    # construction of that legacy shape; confidence callers can still provide
    # the missing fit context explicitly (``confint(..., data=..., family=...)``).
    data: numpy.ndarray = dataclasses.field(
        default_factory=lambda: numpy.empty((0, 0), dtype=float), repr=False
    )
    family: str = 'mean'
    order: tuple = ()
    cp_only: bool = False
    fit_kwargs: Mapping[str, object] = dataclasses.field(
        default_factory=FrozenMapping, repr=False, compare=False
    )
    _copy_arrays: dataclasses.InitVar[bool] = True

    def __post_init__(self, _copy_arrays):
        # Public/manual construction copies every array so later caller
        # mutation cannot change a frozen result. Detector-owned arrays have
        # no external writable aliases and use the private no-copy path to
        # avoid duplicating the full input and every native result buffer.
        # R's S4 result always stores a matrix, so normalize one-dimensional
        # data here too.
        data_value = (
            numpy.empty((0, 0), dtype=float)
            if self.data is None else numpy.asarray(self.data, dtype=float)
        )
        if data_value.ndim == 1:
            data_value = data_value.reshape(-1, 1)
        if data_value.ndim != 2:
            raise ValueError("data must be one- or two-dimensional")
        arrays = {
            'cp_set': numpy.asarray(
                [] if self.cp_set is None else self.cp_set,
                dtype=numpy.int64,
            ),
            'raw_cp_set': numpy.asarray(
                [] if self.raw_cp_set is None else self.raw_cp_set,
                dtype=numpy.int64,
            ),
            'cost_values': numpy.asarray(
                [] if self.cost_values is None else self.cost_values,
                dtype=float,
            ),
            'residuals': numpy.asarray(
                [] if self.residuals is None else self.residuals,
                dtype=float,
            ),
            'thetas': numpy.asarray(
                [] if self.thetas is None else self.thetas,
                dtype=float,
            ),
            'data': numpy.asarray(data_value, dtype=float),
        }
        for name, value in arrays.items():
            if _copy_arrays:
                value = numpy.array(value, copy=True)
            value.setflags(write=False)
            object.__setattr__(self, name, value)
        object.__setattr__(
            self, 'family', 'mean' if self.family is None else str(self.family)
        )
        if self.order is None:
            order_value = ()
        elif isinstance(self.order, str):
            order_value = (self.order,)
        else:
            try:
                order_value = tuple(self.order)
            except TypeError:
                order_value = (self.order,)
        object.__setattr__(self, 'order', order_value)
        object.__setattr__(self, 'cp_only', bool(self.cp_only))
        object.__setattr__(self, 'fit_kwargs', _freeze_mapping(self.fit_kwargs))

    def __reduce__(self):
        # NumPy's pickle restores arrays as writeable and the default slots
        # pickle bypasses __post_init__. Reconstruct through the public
        # initializer so copied/unpickled results retain the read-only
        # contract for arrays and fit options.
        return (
            type(self),
            (
                self.cp_set,
                self.raw_cp_set,
                self.cost_values,
                self.residuals,
                self.thetas,
                self.data,
                self.family,
                self.order,
                self.cp_only,
                self.fit_kwargs,
            ),
        )

    @property
    def details_available(self):
        """Whether costs, residuals, and parameters were computed."""
        return not self.cp_only

    def confint(self, *args, **kwargs):
        """Construct confidence intervals for this result.

        The stored data, family, and order are used automatically.
        """
        from fastcpd.confidence import confint
        return confint(self, *args, **kwargs)

    def summary(self):
        """Return a compact, language-neutral summary of the fit.

        R exposes ``summary.fastcpd`` as a print-oriented method. Returning a
        dictionary keeps the Python API composable while exposing the same
        information to notebooks, logs, and documentation examples.
        """
        try:
            p_response = int(self.fit_kwargs.get('p_response', 0))
        except (TypeError, ValueError, OverflowError):
            p_response = 0
        if p_response > 0:
            n_response = p_response
        elif self.family.lower() in {
            'lm', 'gaussian', 'lasso', 'binomial', 'poisson', 'quantile',
            'ar', 'arma', 'ma', 'arima', 'garch',
        }:
            n_response = 1
        else:
            n_response = self.data.shape[1]
        return {
            'family': self.family,
            'order': self.order,
            'n_obs': int(self.data.shape[0]),
            'n_response': int(n_response),
            'cp_set': self.cp_set.copy(),
            'details_available': self.details_available,
        }

    def plot(self, ax=None, *, column=0, **kwargs):
        """Plot one input column and detected change points.

        Plotting is intentionally an optional Python-layer feature: importing
        fastcpd never requires Matplotlib. The returned Matplotlib axes object
        makes this method convenient in notebooks and mirrors R's plot method.
        """
        try:
            import matplotlib.pyplot as pyplot
        except ImportError as error:  # pragma: no cover - optional dependency
            raise ImportError(
                "CpdResult.plot() requires matplotlib; install it separately"
            ) from error
        if not 0 <= int(column) < self.data.shape[1]:
            raise ValueError("column must refer to an input data column")
        if ax is None:
            _, ax = pyplot.subplots()
        values = self.data[:, int(column)]
        ax.plot(numpy.arange(1, values.size + 1), values, **kwargs)
        for cp in self.cp_set:
            ax.axvline(int(cp), color='grey', linestyle='--', alpha=0.8)
        ax.set_xlabel('Observation')
        ax.set_ylabel('Value')
        return ax


def detect_mean(data, **kwargs):
    """Find change points efficiently in mean change models.

    Args:
        data: Univariate or multivariate data for mean change detection.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A CpdResult. When cp_only=True, detail arrays are empty.
    """
    return detect(data=data, family='mean', **kwargs)


def detect_exponential(data, **kwargs):
    """Find change points efficiently in exponentially distributed data.

    Args:
        data: Univariate data where each observation is exponentially
            distributed; the rate parameter is allowed to change.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A CpdResult. When cp_only=True, detail arrays are empty.
    """
    return detect(data=data, family='exponential', **kwargs)


def detect_variance(data, **kwargs):
    """Find change points efficiently in variance change models.

    Args:
        data: Univariate or multivariate data for variance change detection.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A CpdResult. When cp_only=True, detail arrays are empty.
    """
    return detect(data=data, family='variance', **kwargs)


def detect_meanvariance(data, **kwargs):
    """Find change points efficiently in mean and/or variance change models.

    Args:
        data: Univariate or multivariate data for mean and/or variance change
            detection.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A CpdResult. When cp_only=True, detail arrays are empty.
    """
    return detect(data=data, family='meanvariance', **kwargs)


def detect_var(data, order=0, **kwargs):
    """Find change points efficiently in VAR (vector autoregression) models.

    Args:
        data: Unlagged multivariate time series data, shape (n, q).
        order: Number of lagged predictors per response (p).
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A CpdResult. When cp_only=True, detail arrays are empty.
    """
    return detect(data=data, family='var', order=order, **kwargs)


def detect_lasso(data, **kwargs):
    """Find change points efficiently in LASSO regression models.

    Args:
        data: Data where the first column is the response.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A CpdResult. When cp_only=True, detail arrays are empty.
    """
    return detect(data=data, family='lasso', **kwargs)


def detect_garch(data, order=(0, 0), **kwargs):
    """Find change points in GARCH(p, q) models.

    Args:
        data: Univariate time series, shape (n,) or (n, 1).
        order: Tuple (p, q) — GARCH and ARCH orders.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A CpdResult. When cp_only=True, detail arrays are empty.
    """
    return detect(
        data=_as_univariate_series(data), family='garch', order=order, **kwargs
    )


def detect_lm(data, **kwargs):
    """Find change points in ordinary linear regression models.

    Args:
        data: Array where column 0 is the response and the remaining columns
            are predictors, shape (n, p+1).
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A CpdResult. When cp_only=True, detail arrays are empty.
    """
    return detect(data=data, family='lm', **kwargs)


def detect_binomial(data, **kwargs):
    """Find change points in logistic regression models.

    Args:
        data: Array where column 0 is the binary response and the remaining
            columns are predictors, shape (n, p+1).
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A CpdResult. When cp_only=True, detail arrays are empty.
    """
    return detect(data=data, family='binomial', **kwargs)


def detect_poisson(data, **kwargs):
    """Find change points in Poisson regression models.

    Args:
        data: Array where column 0 is the count response and the remaining
            columns are predictors, shape (n, p+1).
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A CpdResult. When cp_only=True, detail arrays are empty.
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
        A CpdResult. When cp_only=True, detail arrays are empty.
    """
    # Normalize scalar, list, tuple, and NumPy one-element-array spellings
    # through the same validator used by the generic dispatcher.  Comparing
    # a sequence directly with ``0 < order < 1`` raises a TypeError before we
    # can produce the stable ValueError promised by the public API.
    order = _validate_quantile_order(order)
    return detect(data=data, family='quantile', order=(order,), **kwargs)


def detect_arma(data, order=(0, 0), **kwargs):
    """Find change points in ARMA(p, q) models.

    When order[0] == 0 (pure MA), routes to the MA family automatically.

    Args:
        data: Univariate time series, shape (n,) or (n, 1).
        order: Tuple (p, q) — AR and MA orders.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A CpdResult. When cp_only=True, detail arrays are empty.
    """
    return detect(
        data=_as_univariate_series(data), family='arma', order=order, **kwargs
    )


def detect_ar(data, order=0, **kwargs):
    """Find change points in AR(p) models (pure autoregressive).

    Args:
        data: Univariate time series, shape (n,) or (n, 1).
        order: AR order p.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A CpdResult. When cp_only=True, detail arrays are empty.
    """
    return detect(
        data=_as_univariate_series(data), family='ar', order=order, **kwargs
    )


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
        A CpdResult. When cp_only=True, detail arrays are empty.

    For ``d = 0`` this is identical to
    ``arma(data, order=(p, q))``.
    """
    return detect(
        data=_as_univariate_series(data), family='arima', order=order,
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
    data_matrix = _as_2d_data(data)
    transformed = _rank_transform(data_matrix)
    result = detect(data=transformed, family='mean', **kwargs)
    # The public R wrapper reports a mean-family fit after transforming ranks,
    # but bootstrap refits must repeat the rank transform on the original
    # observations.  Keep that distinction in the stored fit metadata.
    rank_fit_kwargs = dict(getattr(result, 'fit_kwargs', {}) or {})
    # Preserve all normalized detector options (including defaults) so a
    # confidence bootstrap can reproduce the original call.  The public R
    # object reports this as a mean-family fit, but the wrapper marker tells
    # Python to repeat the rank transform for every bootstrap sample.
    rank_fit_kwargs['family'] = 'rank'
    rank_fit_kwargs['_wrapper_family'] = 'rank'
    rank_fit_kwargs['_native_family'] = 'mean'
    rank_fit_kwargs['order'] = tuple(result.order)
    return _with_metadata(
        result, data_matrix, family='mean', order=result.order,
        fit_kwargs=rank_fit_kwargs,
    )


def detect_kernel(data, order=(100, 0), random_state=None, **kwargs):
    """Find distributional change points using random Fourier features.

    Args:
        data: Univariate or multivariate data.
        order: Tuple ``(n_features, bandwidth)``. The first component is a
            positive integer number of random features (non-positive values
            use the R-compatible default of 100), and the second is the RBF
            bandwidth. A non-positive bandwidth uses the median heuristic on
            up to 1000 sampled observations. Components after the first two
            are ignored for compatibility with R.
        random_state: Optional NumPy random seed or generator.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A CpdResult. When cp_only=True, detail arrays are empty.
    """
    kwargs.setdefault('cost_adjustment', 'BIC')
    return _detect_kcp_wrapper(
        data, order=order, random_state=random_state, **kwargs
    )


def detect_kcp(data, order=(100, 0), random_state=None, **kwargs):
    """Find distributional change points using the KCP wrapper name."""
    kwargs.setdefault('cost_adjustment', 'BIC')
    return _detect_kcp_wrapper(
        data, order=order, random_state=random_state, **kwargs
    )


def detect_mean_variance(data, **kwargs):
    """Find change points in mean and/or variance change models."""
    return detect(data=data, family='meanvariance', **kwargs)


def detect_linear_regression(data, **kwargs):
    """Find change points in ordinary linear regression models."""
    return detect(data=data, family='lm', **kwargs)


def detect_logistic_regression(data, **kwargs):
    """Find change points in logistic regression models."""
    return detect(data=data, family='binomial', **kwargs)


def detect_poisson_regression(data, **kwargs):
    """Find change points in Poisson regression models."""
    return detect(data=data, family='poisson', **kwargs)


def detect_quantile_regression(data, order=0.5, **kwargs):
    """Find change points in quantile regression models."""
    return detect_quantile(data, order=order, **kwargs)


# R defines these compatibility spellings as direct aliases.  Keep the small
# wrapper definitions above discoverable by the documentation generator, then
# rebind their runtime names so introspection and dispatch have the same
# identity semantics in both languages.
detect_mean_variance = detect_meanvariance  # noqa: F811
detect_linear_regression = detect_lm  # noqa: F811
detect_logistic_regression = detect_binomial  # noqa: F811
detect_poisson_regression = detect_poisson  # noqa: F811
detect_quantile_regression = detect_quantile  # noqa: F811
detect_kcp = detect_kernel  # noqa: F811


mean = detect_mean
exponential = detect_exponential
variance = detect_variance
meanvariance = detect_meanvariance
var = detect_var
mv = detect_meanvariance
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
    cost_adjustment: str = _COST_ADJUSTMENT_UNSET,
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
        cost_adjustment: One of 'BIC', 'MBIC', 'MDL'. Omitting the argument
            uses MBIC; explicitly passing ``None`` selects unadjusted BIC,
            matching R's ``NULL`` behavior.
        family: One of 'mean', 'variance', 'meanvariance', 'exponential',
            'mgaussian' / 'var' (synonym), 'lasso', 'garch', 'gaussian' /
            'lm' (synonym), 'binomial', 'poisson', 'arma', 'ma',
            'quantile', 'arima' (segment-local differencing), or 'kcp'
            (random Fourier feature transform). Use ``detect_rank()`` and
            ``detect_kernel()`` for the wrapper-only R-compatible spellings.
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
        include_mean: Must be False for ARIMA. For other families this
            compatibility option is ignored, as it is in R. The shared ARIMA
            likelihood is zero-mean in both R and Python.
        p: Number of model parameters.  ``None`` (or 0) triggers automatic
            inference from ``family`` and the data dimensions in C++,
            matching the R package's per-family formulas.
        p_response: Number of response columns (mgaussian only).
        variance_estimation: Pre-specified variance/covariance matrix.
            When omitted, mean models use the Rice estimator. Gaussian/LM
            models with a criterion-valued beta and multivariate-regression
            models use R's block-lagged estimator in the shared native layer;
            families that do not consume a covariance use an identity value.
        cp_only: If True, skip segment costs, residuals, and parameter
            estimates. The return type remains ``CpdResult``; its detail
            arrays are empty and ``details_available`` is False.
        vanilla_percentage: Fraction of observations evaluated with pure PELT
            (no gradient update). 1.0 runs full PELT; 0.0 runs full SEN.
        warm_start: If True, use previous segment parameters as initial
            values.
        show_progress: If True, display a tqdm-format progress bar on stderr
            showing PELT timestep progress. Same format as Python tqdm default:
            ``42%|████████████          | 42/100 [00:05<00:07, 8.33it/s]``.
            Implemented in C++; no tqdm package required.

    Returns:
        A ``CpdResult`` with read-only NumPy arrays and fit metadata.
    """
    # Python callers naturally spell the generic API as ``detect(array,
    # family=...)``.  Keep R's formula-first signature for compatibility,
    # while treating a non-string first positional value as the data matrix
    # whenever the explicit ``data=`` argument was omitted.
    if data is None and not isinstance(formula, str):
        data = formula
        formula = 'y ~ . - 1'
    if data is None:
        raise ValueError("data must be provided")
    if multiple_epochs is not None:
        raise NotImplementedError(
            "Custom multiple_epochs schedules are not supported by the "
            "Python binding."
        )

    # Snapshot the user-facing input before constructing lagged designs.  The
    # R result stores the original observations and pads model-lag residuals;
    # confidence refits likewise need this unmodified matrix.
    original_data = _as_2d_data(data)
    data = original_data

    cost_adjustment_missing = cost_adjustment is _COST_ADJUSTMENT_UNSET
    if cost_adjustment_missing:
        cost_adjustment = 'MBIC'
    elif cost_adjustment is None:
        # R's explicit ``cost_adjustment = NULL`` means an unadjusted BIC
        # cost, while omission uses the MBIC default.
        cost_adjustment = 'BIC'
    if not isinstance(cost_adjustment, str):
        raise ValueError(
            "cost_adjustment must be 'BIC', 'MBIC', or 'MDL'"
        )
    if cost_adjustment not in ('BIC', 'MBIC', 'MDL'):
        raise ValueError(
            f"cost_adjustment must be 'BIC', 'MBIC', or 'MDL', "
            f"got {cost_adjustment!r}"
        )

    if isinstance(beta, str):
        if beta not in ('BIC', 'MBIC', 'MDL'):
            raise ValueError("Invalid beta selection criterion provided.")
        beta_value = beta
    else:
        beta_value = _coerce_finite_scalar(beta, "beta")

    trim = _coerce_unit_interval(trim, "trim")
    vanilla_percentage = _coerce_unit_interval(
        vanilla_percentage, "vanilla_percentage"
    )
    segment_count = _coerce_integer(segment_count, "segment_count", minimum=1)
    epsilon = _coerce_finite_scalar(epsilon, "epsilon", positive=True)
    momentum_coef = _coerce_finite_scalar(momentum_coef, "momentum_coef")
    p_response = _coerce_integer(p_response, "p_response", minimum=0)
    if p is not None:
        p = _coerce_integer(p, "p", minimum=0)
    if pruning_coef is not None:
        pruning_coef = _coerce_finite_scalar(pruning_coef, "pruning_coef")
    line_search_array = _validate_option_vector(line_search, "line_search")
    if line_search_array.size and numpy.any(line_search_array <= 0):
        raise ValueError("line_search values must be positive")
    lower_array = _validate_option_vector(lower, "lower")
    upper_array = _validate_option_vector(upper, "upper")
    if (lower_array.size and upper_array.size and
            lower_array.size != upper_array.size):
        raise ValueError("lower and upper must have the same length")
    if lower_array.size and upper_array.size and numpy.any(lower_array > upper_array):
        raise ValueError("lower values must not exceed upper values")

    raw_family = 'custom' if family is None else str(family).lower()
    if (raw_family == 'custom' or cost is not None or
            cost_gradient is not None or cost_hessian is not None):
        raise NotImplementedError(
            "Custom cost callbacks are intentionally R-only; the Python "
            "binding supports the built-in native families only."
        )
    if raw_family not in _SUPPORTED_FAMILIES:
        raise ValueError(
            f"Family '{raw_family}' is not supported by the Python binding. "
            f"Supported families: {sorted(_SUPPORTED_FAMILIES)}."
        )

    public_family = raw_family
    public_order = _public_order(order)
    native_family = raw_family
    native_order = public_order
    native_p = p
    native_p_response = p_response
    native_vanilla = vanilla_percentage
    index_offset = 0

    # KCP is the sole distribution-free family accepted by R's generic
    # dispatcher. Rank and ``kernel`` remain wrapper-only spellings.
    if raw_family == 'kcp':
        # Keep this path in one helper so ``detect_kcp()`` and direct
        # ``detect(family='kcp')`` have exactly the same behavior.
        return _detect_kcp_wrapper(
            original_data, order=public_order, random_state=random_state,
            beta=beta_value,
            cost_adjustment=(
                'BIC' if cost_adjustment_missing else cost_adjustment
            ),
            line_search=line_search_array, lower=lower_array,
            upper=upper_array, pruning_coef=pruning_coef,
            segment_count=segment_count, trim=trim,
            momentum_coef=momentum_coef, epsilon=epsilon,
            p=p, p_response=p_response, cp_only=cp_only,
            warm_start=warm_start, show_progress=show_progress,
            **({'variance_estimation': variance_estimation}
               if variance_estimation is not None else {}),
        )
    elif raw_family == 'var':
        var_order = _validate_var_order(order)
        legacy_columns = p_response * (var_order + 1)
        if p_response > 0 and data.shape[1] == legacy_columns:
            # Compatibility with the pre-1.0 Python ``var`` call shape,
            # where callers supplied [responses, lagged predictors].
            q = p_response
            if q <= 0 or data.shape[1] != q * (var_order + 1):
                raise ValueError("invalid legacy VAR design matrix")
            index_offset = 0
        else:
            data, q = _var_regression_data(data, var_order)
            if p_response not in (0, q):
                raise ValueError(
                    "p_response must match the number of columns in VAR data"
                )
            index_offset = var_order
        native_family = 'mgaussian'
        native_order = (var_order,)
        native_p_response = q
        native_p = var_order * q * q
        native_vanilla = 1.0
    elif raw_family == 'mgaussian':
        # ``mgaussian`` is the advanced spelling for an already constructed
        # [responses, predictors] design.  Unlike the public ``var`` wrapper,
        # it must be told how many leading response columns are present.
        q = p_response
        if q <= 0:
            raise ValueError(
                "mgaussian requires p_response for response/predictor data"
            )
        if q >= data.shape[1]:
            raise ValueError("mgaussian data must include predictors")
        native_p_response = q
        native_p = q * (data.shape[1] - q)
        native_order = (0,)
        native_vanilla = 1.0
    elif raw_family in ('lm', 'gaussian'):
        # R's ``lm`` dispatches multivariate responses to the PELT
        # ``mgaussian`` family.  A scalar response remains the SEN Gaussian
        # route.  ``gaussian`` is retained as an advanced Python spelling.
        if data.shape[1] < 2:
            raise ValueError("regression data must include predictors")
        q = p_response if p_response > 0 else 1
        if q > data.shape[1]:
            raise ValueError("p_response must not exceed data columns")
        if q > 1:
            if data.shape[1] <= q:
                raise ValueError("multivariate regression requires predictors")
            native_family = 'mgaussian'
            native_p_response = q
            native_p = q * (data.shape[1] - q)
            native_vanilla = 1.0
        else:
            native_family = 'gaussian'
            native_p_response = 1
            if native_p is None:
                native_p = data.shape[1] - 1
    elif raw_family == 'ar':
        ar_order = _validate_ar_order(order)
        data = _make_ar_design(data, ar_order)
        native_family = 'gaussian'
        native_order = (ar_order,)
        native_p_response = 1
        native_p = ar_order
        index_offset = ar_order
    elif raw_family == 'arima':
        if not isinstance(include_mean, (bool, numpy.bool_)):
            raise ValueError("include_mean must be a single logical value")
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
        native_p = p_ar + q_ma + 1
        if d_int == 0:
            native_order = (p_ar, q_ma)
            if p_ar == 0:
                native_family = 'ma'
                native_p = q_ma + 1
            elif q_ma == 0:
                data = _make_ar_design(data, p_ar)
                native_family = 'gaussian'
                native_order = (p_ar,)
                index_offset = p_ar
                native_p = p_ar
            else:
                native_family = 'arma'
        else:
            native_family = 'arima'
            native_order = (p_ar, d_int, q_ma)
            native_vanilla = 1.0
    elif raw_family == 'arma':
        arma_order = _validate_arma_order(order)
        p_ar, q_ma = arma_order
        if data.shape[1] != 1:
            raise ValueError("ARMA data must be univariate")
        native_order = arma_order
        if p_ar == 0:
            native_family = 'ma'
            native_p = q_ma + 1
        elif q_ma == 0:
            data = _make_ar_design(data, p_ar)
            native_family = 'gaussian'
            native_order = (p_ar,)
            native_p = p_ar
            index_offset = p_ar
        else:
            native_family = 'arma'
            native_p = p_ar + q_ma + 1
    elif raw_family == 'ma':
        native_order = _validate_ma_order(order)
        native_p = native_order[1] + 1
    elif raw_family == 'garch':
        if data.shape[1] != 1:
            raise ValueError("GARCH data must be univariate")
        native_order = _validate_garch_order(order)
        native_p = sum(native_order) + 1
        native_vanilla = 1.0
    elif raw_family == 'quantile':
        if data.shape[1] < 2:
            raise ValueError("quantile data must include predictors")
        native_order = (_validate_quantile_order(order),)
    elif raw_family == 'mean':
        native_order = (0, 0, 0)
        native_p = data.shape[1]
        native_vanilla = 1.0
    elif raw_family == 'variance':
        native_order = (0, 0, 0)
        native_p = data.shape[1] ** 2
        native_vanilla = 1.0
    elif raw_family == 'meanvariance':
        native_order = (0, 0, 0)
        native_p = data.shape[1] + data.shape[1] ** 2
        native_vanilla = 1.0
    elif raw_family == 'exponential':
        native_order = (0, 0, 0)
        native_p = data.shape[1]
        native_vanilla = 1.0

    if raw_family in ('lasso', 'binomial', 'poisson') and data.shape[1] < 2:
        raise ValueError(f"{raw_family} data must include predictors")

    # R's character beta path multiplies the scalar Gaussian penalty by
    # ``c(sigma_)``.  A supplied covariance with more than one element then
    # fails at the Rcpp scalar conversion.  Reject the same ambiguous shape
    # before crossing the Python/native boundary; numeric beta remains
    # permissive just as in R because no scalar rescaling is attempted.
    if (isinstance(beta_value, str) and native_family == 'gaussian' and
            variance_estimation is not None):
        supplied_variance = _coerce_variance_matrix(variance_estimation)
        if supplied_variance.shape != (1, 1):
            raise ValueError(
                "variance_estimation for Gaussian/lm fits must be scalar "
                "when beta is a character criterion"
            )

    # ``include_mean`` is an ellipsis-only compatibility option in R and is
    # intentionally ignored for every non-ARIMA family.
    native_data = numpy.ascontiguousarray(data, dtype=float)
    if (variance_estimation is None and
            native_family in ('mean', 'gaussian', 'mgaussian')):
        # Let the shared C++ core perform the R-compatible automatic
        # estimator.  Passing an empty 0×0 matrix is the pybind sentinel for
        # omission; computing the estimate here would bypass the common
        # block-lagged variance.lm()/nearest-PD implementation.
        variance_value = numpy.empty((0, 0), dtype=float)
    else:
        variance_value = _estimate_variance(
            native_data, native_family, native_p_response,
            variance_estimation
        )
        variance_value = _coerce_variance_matrix(variance_value)

    # p is forced by R for fixed-dimensional PELT/time-series routes.  For
    # ordinary regression/SEN families an explicit p remains meaningful.
    p_int = int(native_p) if native_p is not None else 0
    pruning_float = (
        float('nan') if pruning_coef is None else float(pruning_coef)
    )

    result = fastcpd_impl(
        beta_value,
        cost_adjustment,
        bool(cp_only),
        native_data,
        epsilon,
        native_family,
        line_search_array,
        lower_array,
        momentum_coef,
        numpy.asarray(native_order, dtype=float),
        p_int,
        int(native_p_response),
        pruning_float,
        segment_count,
        trim,
        upper_array,
        native_vanilla,
        variance_value,
        bool(warm_start),
        bool(show_progress),
    )

    residuals = numpy.asarray(result['residuals'], dtype=float)
    if not cp_only and index_offset > 0:
        if residuals.ndim == 1:
            residuals = residuals.reshape(-1, 1)
        if residuals.ndim == 2:
            residuals = numpy.vstack([
                numpy.full((index_offset, residuals.shape[1]), numpy.nan),
                residuals,
            ])

    fit_kwargs = _build_fit_kwargs(
        family=public_family,
        order=public_order,
        beta=beta_value,
        cost_adjustment=cost_adjustment,
        line_search=line_search_array,
        lower=lower_array if lower_array.size else None,
        upper=upper_array if upper_array.size else None,
        pruning_coef=pruning_coef,
        segment_count=segment_count,
        trim=trim,
        momentum_coef=momentum_coef,
        epsilon=epsilon,
        p=p,
        p_response=p_response,
        variance_estimation=variance_estimation,
        vanilla_percentage=native_vanilla,
        warm_start=warm_start,
        show_progress=show_progress,
        random_state=random_state if raw_family == 'kcp' else None,
    )

    cp_set = numpy.asarray(result['cp_set'], dtype=numpy.int64)
    raw_cp_set = numpy.asarray(result['raw_cp_set'], dtype=numpy.int64)
    if index_offset:
        cp_set = cp_set + index_offset
        raw_cp_set = raw_cp_set + index_offset

    return CpdResult(
        cp_set=cp_set,
        raw_cp_set=raw_cp_set,
        cost_values=result['cost_values'],
        residuals=residuals,
        thetas=result['thetas'],
        data=original_data,
        family=public_family,
        order=public_order,
        cp_only=cp_only,
        fit_kwargs=fit_kwargs,
        _copy_arrays=False,
    )


def _public_order(order):
    """Return model order metadata as a tuple without changing its values."""
    if isinstance(order, numpy.ndarray) and order.ndim == 0:
        return (order.item(),)
    if hasattr(order, '__len__') and not isinstance(order, str):
        try:
            return tuple(order)
        except TypeError:
            pass
    return (order,)


def _optional_vector(value):
    """Convert an optional vector argument to a contiguous NumPy buffer."""
    if value is None:
        return numpy.empty(0, dtype=float)
    array = numpy.atleast_1d(numpy.asarray(value, dtype=float))
    if array.ndim != 1:
        raise ValueError("option values must be one-dimensional")
    return numpy.ascontiguousarray(array)


def _coerce_integer(value, name, *, minimum=None):
    """Coerce a scalar to an integer while rejecting lossy values."""
    try:
        value_float = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be an integer") from error
    if not numpy.isfinite(value_float) or value_float != numpy.floor(value_float):
        raise ValueError(f"{name} must be an integer")
    value_int = int(value_float)
    if minimum is not None and value_int < minimum:
        if minimum == 1:
            raise ValueError(f"{name} must be a positive integer")
        raise ValueError(f"{name} must be at least {minimum}")
    return value_int


def _coerce_unit_interval(value, name):
    try:
        value_float = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be a finite value in [0, 1]") from error
    if not numpy.isfinite(value_float) or not 0 <= value_float <= 1:
        raise ValueError(f"{name} must be a finite value in [0, 1]")
    return value_float


def _coerce_finite_scalar(value, name, *, positive=False):
    try:
        value_float = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError(f"{name} must be a finite numeric value") from error
    if not numpy.isfinite(value_float) or (positive and value_float <= 0):
        suffix = " greater than zero" if positive else ""
        raise ValueError(f"{name} must be a finite numeric value{suffix}")
    return value_float


def _validate_option_vector(value, name):
    """Validate lower/upper/line-search vectors before crossing pybind."""
    array = _optional_vector(value)
    if array.size and numpy.any(numpy.isnan(array)):
        raise ValueError(f"{name} must not contain NaN")
    return array


def _make_ar_design(data, ar_order):
    """Construct R's lagged [response, predictors] AR design matrix."""
    if data.shape[1] != 1:
        raise ValueError("AR data must be univariate")
    n_rows = data.shape[0]
    if n_rows <= ar_order:
        raise ValueError("AR order must be smaller than the number of rows")
    response = data[ar_order:, 0:1]
    lags = numpy.column_stack([
        data[ar_order - lag:n_rows - lag, 0]
        for lag in range(1, ar_order + 1)
    ])
    return numpy.column_stack([response, lags])


def _coerce_variance_matrix(value):
    """Return a two-dimensional covariance buffer for pybind."""
    matrix = numpy.asarray(value, dtype=float)
    if matrix.ndim == 0:
        matrix = matrix.reshape(1, 1)
    if matrix.ndim != 2:
        raise ValueError("variance_estimation must be a scalar or 2-D matrix")
    if matrix.shape[0] != matrix.shape[1] or matrix.shape[0] == 0:
        raise ValueError("variance_estimation must be a non-empty square matrix")
    if not numpy.all(numpy.isfinite(matrix)):
        raise ValueError("variance_estimation must contain only finite values")
    return numpy.ascontiguousarray(matrix, dtype=float)


def _integer_order_values(order, expected_length, family):
    """Return a fixed-length tuple of non-negative integer order values."""
    if not hasattr(order, '__len__') or isinstance(order, str):
        raise ValueError(
            f"{family} order must contain {expected_length} integers"
        )
    try:
        values = list(order)
    except TypeError as error:
        raise ValueError(
            f"{family} order must contain {expected_length} integers"
        ) from error
    if len(values) != expected_length:
        raise ValueError(
            f"{family} order must contain {expected_length} integers"
        )
    integers = []
    for value in values:
        try:
            value_float = float(value)
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError(
                f"{family} order must contain non-negative integers"
            ) from error
        value_int = int(value_float) if numpy.isfinite(value_float) else -1
        if (not numpy.isfinite(value_float) or value_int < 0 or
                value_float != value_int):
            raise ValueError(
                f"{family} order must contain non-negative integers"
            )
        integers.append(value_int)
    return tuple(integers)


def _validate_ar_order(order):
    """Return a validated positive AR order.

    R accepts both ``p`` and the three-component ARIMA-style spelling
    ``c(p, 0, 0)`` for the AR wrapper.  The latter was historically rejected
    by Python even though it describes the same model.
    """
    if hasattr(order, '__len__') and not isinstance(order, str):
        try:
            values = list(order)
        except TypeError as error:
            raise ValueError("AR order must be a positive integer") from error
        if len(values) == 3:
            values = _integer_order_values(values, 3, "AR")
            if values[0] <= 0 or values[1:] != (0, 0):
                raise ValueError(
                    "AR order must be p or (p, 0, 0), with p positive"
                )
            return values[0]
        if len(values) != 1:
            raise ValueError("AR order must be a positive integer")
        order = values[0]
    return _coerce_integer(order, "AR order", minimum=1)


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
        try:
            order_values = list(order)
        except TypeError as error:
            raise ValueError("VAR order must be a positive integer") from error
        if len(order_values) != 1:
            raise ValueError("VAR order must be a positive integer")
        order = order_values[0]
    return _coerce_integer(order, "VAR order", minimum=1)


def _validate_garch_order(order):
    values = _integer_order_values(order, 2, "GARCH")
    if values == (0, 0):
        raise ValueError("GARCH order must contain at least one non-zero value")
    return values


def _validate_ma_order(order):
    """Normalize MA order to the native ``(0, q)`` spelling."""
    if not hasattr(order, '__len__') or isinstance(order, str):
        q = _coerce_integer(order, "MA order", minimum=1)
        return (0, q)
    values = _integer_order_values(order, 2, "MA")
    if values[0] != 0 or values[1] <= 0:
        raise ValueError("MA order must be (0, q), with q positive")
    return values


def _validate_quantile_order(order):
    if hasattr(order, '__len__') and not isinstance(order, str):
        try:
            values = list(order)
        except TypeError as error:
            raise ValueError("quantile order must contain one level") from error
        if len(values) != 1:
            raise ValueError("quantile order must contain one level")
        order = values[0]
    value = _coerce_finite_scalar(order, "quantile order")
    if not 0 < value < 1:
        raise ValueError(f"order must be in (0, 1), got {order!r}")
    return value


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


def _validate_kernel_order(order):
    """Normalize the R-compatible KCP order before allocating features.

    R consumes at most the first two entries of ``order``.  An omitted or
    empty order uses 100 random features and the median bandwidth heuristic;
    non-positive feature counts likewise fall back to 100, while a
    non-positive bandwidth selects the heuristic.  Python rejects
    non-finite and positive fractional values early so malformed input cannot
    reach a zero-width or non-finite random-feature matrix.
    """
    if order is None:
        values = []
    elif isinstance(order, numpy.ndarray) and order.ndim == 0:
        values = [order.item()]
    elif isinstance(order, str):
        values = [order]
    else:
        try:
            values = list(order)
        except TypeError:
            values = [order]

    # R's implementation only reads order[1] and order[2].  Keep that
    # behavior rather than making otherwise-valid legacy calls fail because
    # they carry metadata after the bandwidth.
    values = values[:2]
    if len(values) == 1 and values[0] is None:
        values = []

    feature_count = 100
    if values:
        try:
            feature_value = float(values[0])
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError(
                "KCP order[0] must be a finite numeric feature count"
            ) from error
        if not numpy.isfinite(feature_value):
            raise ValueError(
                "KCP order[0] must be a finite numeric feature count"
            )
        if feature_value > 0:
            if feature_value != numpy.floor(feature_value):
                raise ValueError(
                    "KCP order[0] must be a positive integer when positive"
                )
            feature_count = int(feature_value)
            # R's integer vectors cannot represent a larger feature count;
            # reject it before an accidentally enormous allocation in Python.
            if feature_count > numpy.iinfo(numpy.int32).max:
                raise ValueError("KCP feature count is too large")

    bandwidth = 0.0
    if len(values) >= 2 and values[1] is not None:
        try:
            bandwidth = float(values[1])
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError(
                "KCP order[1] must be a finite numeric bandwidth"
            ) from error
        if not numpy.isfinite(bandwidth):
            raise ValueError(
                "KCP order[1] must be a finite numeric bandwidth"
            )
    return feature_count, bandwidth


def _kernel_transform(data, order=(100, 0), random_state=None):
    data_matrix = numpy.asarray(data, dtype=float)
    if data_matrix.ndim == 1:
        data_matrix = data_matrix.reshape(-1, 1)
    feature_count, bandwidth = _validate_kernel_order(order)
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
        return _coerce_variance_matrix(variance_estimation)
    if family == 'mean':
        if data.shape[0] < 2:
            return numpy.eye(data.shape[1])
        return _variance_estimation.estimate_variance_mean(data)
    if family == 'mgaussian':
        q = p_response if p_response > 0 else data.shape[1]
        if q > data.shape[1]:
            raise ValueError("p_response must not exceed data columns")
        if data.shape[1] > q:
            # Match R's block-lagged variance.lm() estimator.  It is more
            # expensive than a single OLS Rice estimate, but is the contract
            # used for VAR and multivariate lm models in the R package.
            estimate = _variance_estimation.estimate_variance_linear_regression(
                data, d=q
            )
            if numpy.ndim(estimate) == 0 or numpy.any(~numpy.isfinite(estimate)):
                # Singular short blocks are common in tiny examples.  R's
                # rcond guard ultimately falls back to a tiny diagonal matrix;
                # use identity here so pybind receives a valid covariance.
                return numpy.eye(q)
            return estimate
        if data.shape[0] < 2:
            return numpy.eye(q)
        return _variance_estimation.estimate_variance_mean(
            data[:, :q]
        )
    if family == 'gaussian':
        # Scalar lm/AR fits use the same variance.lm() estimator as R.  A
        # robust fallback keeps short or singular designs usable.
        try:
            estimate = _variance_estimation.estimate_variance_linear_regression(
                data, d=1
            )
            if numpy.ndim(estimate) == 0 and numpy.isfinite(estimate):
                return numpy.asarray([[float(estimate)]])
        except (ValueError, numpy.linalg.LinAlgError):
            pass
        return numpy.eye(1)
    # All other families: variance_estimate is not used by the C++ cost
    # function, so a 1×1 identity placeholder is sufficient.
    return numpy.eye(1)


# R-compatible convenience spellings.  Keeping these as direct aliases (and
# exporting them from ``fastcpd``) makes examples portable without introducing
# a second dispatcher or changing the canonical ``detect_*`` names.
fastcpd = detect
fastcpd_mean = detect_mean
fastcpd_variance = detect_variance
fastcpd_meanvariance = detect_meanvariance
fastcpd_mean_variance = detect_mean_variance
fastcpd_mv = detect_meanvariance
fastcpd_exponential = detect_exponential
fastcpd_var = detect_var
fastcpd_lasso = detect_lasso
fastcpd_garch = detect_garch
fastcpd_lm = detect_lm
fastcpd_linear_regression = detect_linear_regression
fastcpd_binomial = detect_binomial
fastcpd_logistic_regression = detect_logistic_regression
fastcpd_poisson = detect_poisson
fastcpd_poisson_regression = detect_poisson_regression
fastcpd_quantile = detect_quantile
fastcpd_quantile_regression = detect_quantile_regression
fastcpd_arma = detect_arma
fastcpd_ar = detect_ar
fastcpd_arima = detect_arima
fastcpd_rank = detect_rank
fastcpd_kernel = detect_kernel
fastcpd_kcp = detect_kcp
