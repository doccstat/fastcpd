"""
Variance estimation for change point detection models.
"""

import collections
import numpy

VarianceArmaResult = collections.namedtuple(
    'VarianceArmaResult',
    ['table', 'sigma2_aic', 'sigma2_bic'],
)


def estimate_variance(data, family='mean', **kwargs):
    """Estimate variance using the named model family."""
    family = family.lower().replace('-', '_').replace('.', '_')
    if family == 'mean':
        return estimate_variance_mean(data, **kwargs)
    if family == 'median':
        return estimate_variance_median(data, **kwargs)
    if family in ('linear_regression', 'lm'):
        return estimate_variance_linear_regression(data, **kwargs)
    if family == 'arma':
        return estimate_variance_arma(data, **kwargs)
    raise ValueError(
        "family must be one of 'mean', 'median', 'linear_regression', "
        "'lm', or 'arma'"
    )


def _as_numeric_array(data, *, ndim=None, name="data"):
    """Return a finite floating-point array with a predictable shape."""
    values = numpy.asarray(data, dtype=float)
    if ndim is not None and values.ndim != ndim:
        raise ValueError(f"{name} must be a {ndim}-D numeric array")
    if values.size == 0:
        raise ValueError(f"{name} must be non-empty")
    if not numpy.all(numpy.isfinite(values)):
        raise ValueError(f"{name} must contain only finite numeric values")
    return values


def estimate_variance_mean(data):
    """
    Variance estimation for mean change models (Rice estimator).

    data : array-like, shape (n, p)
      Each row is a p-vector observation.

    Returns
    -------
    ndarray, shape (p, p)
      Estimated variance-covariance matrix.
    """
    data_matrix = _as_numeric_array(data)
    if data_matrix.ndim == 0 or data_matrix.ndim > 2:
        raise ValueError("data must be one- or two-dimensional")
    if data_matrix.ndim == 1:
        data_matrix = data_matrix.reshape(-1, 1)
    if data_matrix.shape[0] < 2:
        return numpy.full((data_matrix.shape[1], data_matrix.shape[1]), numpy.nan)
    diffs = data_matrix[1:] - data_matrix[:-1]
    return numpy.mean(diffs[:, :, None] * diffs[:, None, :], axis=0) / 2


def estimate_variance_median(data):
    """
    Variance estimation for median change models (Rice estimator).

    data : array-like, shape (n,)
      Univariate series.

    Returns
    -------
    float
      Estimated variance.
    """
    data_flat = _as_numeric_array(data).ravel(order='F')
    if data_flat.size < 2:
        return float('nan')
    return 2 * (2 * numpy.mean(numpy.abs(numpy.diff(data_flat))) / 3) ** 2


def estimate_variance_linear_regression(
    data,
    d=1,
    block_size=None,
    outlier_iqr=numpy.inf,
):
    """Estimate residual variance for piecewise linear regression models."""
    data_matrix = _as_numeric_array(data, ndim=2)
    # Validate ``d`` before deriving the default block size.  Computing
    # ``ncol(data) - d + 1`` first used to leak TypeError/OverflowError for
    # malformed values (for example ``d=None`` or ``d=inf``), and could derive
    # a non-integral default for values such as ``d=1.5``.  R validates the
    # response dimension as part of the estimator contract, so Python should
    # report one stable ValueError before doing any arithmetic with it.
    try:
        d_float = float(d)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("d must be an integer") from error
    if (not numpy.isfinite(d_float) or d_float != numpy.floor(d_float)):
        raise ValueError("d must be an integer")
    d = int(d_float)
    if d <= 0 or d > data_matrix.shape[1]:
        raise ValueError("d must be between 1 and the number of columns")

    if block_size is None:
        block_size = data_matrix.shape[1] - d + 1
    try:
        block_size_float = float(block_size)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("block_size must be an integer") from error
    if (not numpy.isfinite(block_size_float) or
            block_size_float != numpy.floor(block_size_float)):
        raise ValueError("block_size must be an integer")
    block_size = int(block_size_float)
    n_rows = data_matrix.shape[0]
    if block_size <= 0 or block_size >= n_rows:
        raise ValueError("block_size must be positive and smaller than nrow(data)")

    estimators = numpy.full((n_rows - block_size, d, d), numpy.nan)
    for i in range(n_rows - block_size):
        block = slice(i, i + block_size)
        block_lagged = slice(i + 1, i + block_size + 1)
        y_block = data_matrix[block, :d]
        x_block = data_matrix[block, d:]
        y_block_lagged = data_matrix[block_lagged, :d]
        x_block_lagged = data_matrix[block_lagged, d:]
        try:
            x_t_x = x_block.T @ x_block
            x_t_x_lagged = x_block_lagged.T @ x_block_lagged
            # ``base::solve.default`` rejects a matrix when LAPACK's
            # reciprocal one-norm condition estimate is below its default
            # tolerance (machine epsilon).  ``numpy.linalg.solve`` often
            # proceeds for the same nearly singular matrix and can produce a
            # large, negative variance estimate.  Apply the corresponding
            # guard explicitly so failed blocks are omitted in both
            # languages.
            if (_reciprocal_condition_1norm(x_t_x) < numpy.finfo(float).eps or
                    _reciprocal_condition_1norm(x_t_x_lagged) <
                    numpy.finfo(float).eps):
                raise numpy.linalg.LinAlgError("computationally singular block")
            block_slope = numpy.linalg.solve(x_t_x, x_block.T @ y_block)
            block_lagged_slope = numpy.linalg.solve(
                x_t_x_lagged,
                x_block_lagged.T @ y_block_lagged,
            )
            x_t_x_inv = numpy.linalg.inv(x_t_x)
            x_t_x_inv_lagged = numpy.linalg.inv(x_t_x_lagged)
            cross_term_x = x_block[1:, :].T @ x_block_lagged[:-1, :]
            cross_term = x_t_x_inv @ x_t_x_inv_lagged @ cross_term_x
            slope_delta = block_slope - block_lagged_slope
            delta_numerator = slope_delta.T @ slope_delta
            delta_denominator = numpy.zeros((d, d))
            for j in range(d):
                for k in range(d):
                    if j != k:
                        delta = block_slope[:, j] - block_lagged_slope[:, k]
                        delta_denominator[j, k] += delta.T @ delta
            delta_denominator = (
                delta_denominator +
                numpy.trace(x_t_x_inv + x_t_x_inv_lagged - 2 * cross_term)
            )
            estimators[i, :, :] = delta_numerator / delta_denominator
        except numpy.linalg.LinAlgError:
            continue

    if d == 1:
        values = estimators.ravel()
        values = values[~numpy.isnan(values)]
        if values.size == 0:
            return numpy.nan
        with numpy.errstate(invalid='ignore'):
            q75 = numpy.quantile(values, 0.75)
            iqr = numpy.quantile(values, 0.75) - numpy.quantile(values, 0.25)
            outlier_threshold = q75 + outlier_iqr * iqr
        # R evaluates the threshold even for its default ``Inf`` multiplier.
        # In particular, a zero IQR produces ``Inf * 0 = NaN`` and therefore
        # an empty retained set whose mean is NaN, rather than zero.
        retained = values[values < outlier_threshold]
        return numpy.mean(retained) if retained.size else numpy.nan
    # R's ``colMeans(..., na.rm = TRUE)`` on the (n, d, d) array averages
    # over the block dimension and preserves the d-by-d layout.  Avoid
    # ``numpy.nanmean`` here: it emits ``RuntimeWarning: Mean of empty slice``
    # when every solve for a cell failed, while R quietly returns ``NaN``.
    valid = ~numpy.isnan(estimators)
    counts = numpy.sum(valid, axis=0)
    with numpy.errstate(invalid='ignore', divide='ignore'):
        sums = numpy.sum(numpy.where(valid, estimators, 0.0), axis=0)
        return numpy.divide(
            sums,
            counts,
            out=numpy.full((d, d), numpy.nan),
            where=counts > 0,
        )


def _reciprocal_condition_1norm(matrix):
    """Approximate R/LAPACK's reciprocal one-norm condition estimate."""
    try:
        condition = numpy.linalg.cond(matrix, p=1)
    except numpy.linalg.LinAlgError:
        return 0.0
    if not numpy.isfinite(condition) or condition <= 0:
        return 0.0
    return float(1.0 / condition)


def estimate_variance_arma(data, p, q, max_order=None):
    """Estimate innovation variance for ARMA models via AR approximations."""
    data_flat = _as_numeric_array(data).ravel(order='F')
    try:
        p_float, q_float = float(p), float(q)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("p and q must be non-negative integers") from error
    p_int, q_int = int(p_float), int(q_float)
    if (p_float != p_int or q_float != q_int or p_int < 0 or q_int < 0):
        raise ValueError("p and q must be non-negative integers")
    if max_order is None:
        max_order = p_int * q_int
    try:
        max_order_float = float(max_order)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("max_order must be a positive integer") from error
    max_order = int(max_order_float)
    if max_order_float != max_order:
        raise ValueError("max_order must be a positive integer")
    if max_order <= 0:
        raise ValueError("max_order must be positive")
    if data_flat.size <= max_order + 1:
        raise ValueError("data must contain more observations than max_order")

    rows = []
    for order in range(1, max_order + 1):
        y = data_flat[order:]
        x = numpy.column_stack([
            data_flat[order - lag - 1:len(data_flat) - lag - 1]
            for lag in range(order)
        ])
        sigma2 = estimate_variance_linear_regression(
            numpy.column_stack([y, x])
        )
        # Keep the R result's NA behavior for degenerate blocks while making
        # the information criteria explicitly non-finite in that case.
        if not numpy.isfinite(sigma2) or sigma2 <= 0:
            aic = bic = float('nan')
        else:
            aic = numpy.log(sigma2) + 2 * order / data_flat.size
            bic = (
                numpy.log(sigma2) +
                order * numpy.log(data_flat.size) / data_flat.size
            )
        rows.append({
            'model': f'AR({order})',
            'sigma2': sigma2,
            'AIC': aic,
            'BIC': bic,
        })

    aic_values = numpy.asarray([row['AIC'] for row in rows])
    bic_values = numpy.asarray([row['BIC'] for row in rows])
    if not numpy.any(numpy.isfinite(aic_values)):
        sigma2_aic = float('nan')
    else:
        sigma2_aic = rows[int(numpy.nanargmin(aic_values))]['sigma2']
    if not numpy.any(numpy.isfinite(bic_values)):
        sigma2_bic = float('nan')
    else:
        sigma2_bic = rows[int(numpy.nanargmin(bic_values))]['sigma2']
    return VarianceArmaResult(
        table=rows,
        sigma2_aic=sigma2_aic,
        sigma2_bic=sigma2_bic,
    )


def estimate_variance_lm(data, *args, **kwargs):
    """Estimate residual variance for linear regression models."""
    return estimate_variance_linear_regression(data, *args, **kwargs)


# R defines ``estimate_variance_lm`` as a direct alias of
# ``estimate_variance_linear_regression``.  Keep the wrapper source above
# discoverable to the documentation generator, but expose the same callable
# identity at runtime for portable alias checks and dispatch.
estimate_variance_lm = estimate_variance_linear_regression  # noqa: F811


# R-compatible descriptive aliases.  The dotted R spellings cannot be
# identifiers in Python, so the underscore forms are the canonical names.
variance_mean = estimate_variance_mean
variance_median = estimate_variance_median
variance_lm = estimate_variance_lm
variance_linear_regression = estimate_variance_linear_regression
variance_arma = estimate_variance_arma


mean = estimate_variance_mean
median = estimate_variance_median
lm = estimate_variance_lm
arma = estimate_variance_arma


__all__ = [
    'VarianceArmaResult',
    'arma',
    'estimate_variance',
    'estimate_variance_arma',
    'estimate_variance_linear_regression',
    'estimate_variance_lm',
    'estimate_variance_mean',
    'estimate_variance_median',
    'lm',
    'mean',
    'median',
    'variance_arma',
    'variance_linear_regression',
    'variance_lm',
    'variance_mean',
    'variance_median',
]
