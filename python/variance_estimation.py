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
    data_matrix = numpy.asarray(data)
    if data_matrix.ndim == 1:
        data_matrix = data_matrix.reshape(-1, 1)
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
    data_flat = numpy.asarray(data).ravel()
    return 2 * (2 * numpy.mean(numpy.abs(numpy.diff(data_flat))) / 3) ** 2


def estimate_variance_linear_regression(
    data,
    d=1,
    block_size=None,
    outlier_iqr=numpy.inf,
):
    """Estimate residual variance for piecewise linear regression models."""
    data_matrix = numpy.asarray(data, dtype=float)
    if data_matrix.ndim != 2:
        raise ValueError("data must be a 2-D array with responses first")
    if block_size is None:
        block_size = data_matrix.shape[1] - d + 1
    block_size = int(block_size)
    d = int(d)
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
        if numpy.isfinite(outlier_iqr):
            q75 = numpy.quantile(values, 0.75)
            iqr = numpy.quantile(values, 0.75) - numpy.quantile(values, 0.25)
            values = values[values < q75 + outlier_iqr * iqr]
        return numpy.mean(values)
    return numpy.nanmean(estimators, axis=0)


def estimate_variance_arma(data, p, q, max_order=None):
    """Estimate innovation variance for ARMA models via AR approximations."""
    data_flat = numpy.asarray(data, dtype=float).ravel()
    if max_order is None:
        max_order = int(p) * int(q)
    max_order = int(max_order)
    if max_order <= 0:
        raise ValueError("max_order must be positive")

    rows = []
    for order in range(1, max_order + 1):
        y = data_flat[order:]
        x = numpy.column_stack([
            data_flat[order - lag - 1:len(data_flat) - lag - 1]
            for lag in range(order)
        ])
        sigma2 = estimate_variance_linear_regression(numpy.column_stack([y, x]))
        rows.append({
            'model': f'AR({order})',
            'sigma2': sigma2,
            'AIC': numpy.log(sigma2) + 2 * order / len(data_flat),
            'BIC': (
                numpy.log(sigma2) +
                order * numpy.log(len(data_flat)) / len(data_flat)
            ),
        })

    aic_values = numpy.asarray([row['AIC'] for row in rows])
    bic_values = numpy.asarray([row['BIC'] for row in rows])
    return VarianceArmaResult(
        table=rows,
        sigma2_aic=rows[int(numpy.nanargmin(aic_values))]['sigma2'],
        sigma2_bic=rows[int(numpy.nanargmin(bic_values))]['sigma2'],
    )


def estimate_variance_lm(data, *args, **kwargs):
    """Estimate residual variance for linear regression models."""
    return estimate_variance_linear_regression(data, *args, **kwargs)


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
]
