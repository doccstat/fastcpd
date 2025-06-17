"""
Perform change point detection using fastcpd.
"""

# import pandas as pd
import numpy
# from patsy import dmatrices

import fastcpd.variance_estimation
from math import log, log2
from fastcpd.interface import fastcpd_impl


def mean(data, **kwargs):
    """Find change points efficiently in mean change models.

    Args:
        data: Univariate or multivariate data for mean change detection.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A ``FastCPDResult`` object.
    """
    return detect(data=data, family='mean', **kwargs)


def variance(data, **kwargs):
    """Find change points efficiently in variance change models.

    Args:
        data: Univariate or multivariate data for variance change detection.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A ``FastCPDResult`` object.
    """
    return detect(data=data, family='variance', **kwargs)


def meanvariance(data, **kwargs):
    """Find change points efficiently in mean and/or variance change models.

    Args:
        data: Univariate or multivariate data for mean and/or variance change
            detection.
        **kwargs: Additional arguments passed to ``detect()``.

    Returns:
        A ``FastCPDResult`` object.
    """
    return detect(data=data, family='meanvariance', **kwargs)


def detect(
    formula: str = 'y ~ . - 1',
    data: numpy.ndarray = None,
    beta: object = 'MBIC',
    cost_adjustment: str = 'MBIC',
    family: str = None,
    cost=None,
    cost_gradient=None,
    cost_hessian=None,
    line_search=(1,),
    lower=None,
    upper=None,
    pruning_coef: float = None,
    segment_count: int = 10,
    trim: float = 0.05,
    momentum_coef: float = 0,
    multiple_epochs=lambda x: 0,
    epsilon: float = 1e-10,
    order=(0, 0, 0),
    p: int = None,
    variance_estimation=None,
    cp_only: bool = False,
    vanilla_percentage: float = 0,
    warm_start: bool = False,
    **kwargs
):
    r"""Find change points efficiently.

    Args:
        formula: A formula string specifying the model to be fitted. The
            optional response variable should be on the LHS, covariates on the
            RHS. Intercept should be removed by appending '- 1'. By default,
            an intercept column is added as in R's lm().
        data: A NumPy array of shape (T, d) containing the data to be
            segmented. Each row is a data point $z_t$ in $\mathbb{R}^d$.
        beta: Penalty criterion for the number of change points. Can be one of
            'BIC', 'MBIC', 'MDL', or a float value.
        cost_adjustment: Cost adjustment criterion modifying the cost function.
            Can be one of 'BIC', 'MBIC', 'MDL', or None.
        family: Family of change point model. One of: 'mean', 'variance',
            'meanvariance', 'lm', 'binomial', 'poisson', 'lasso', 'ar',
            'arma', 'arima', 'garch', 'var', 'custom'. If None, it is
            treated as 'custom'.
        cost: Custom cost function, e.g., ``cost(data)`` or
            ``cost(data, theta)``.
        cost_gradient: Gradient of custom cost, e.g.,
            ``cost_gradient(data, theta)``.
        cost_hessian: Hessian of custom cost, e.g.,
            ``cost_hessian(data, theta)``.
        line_search: Values for line search step sizes.
        lower: Lower bound for parameters after each update.
        upper: Upper bound for parameters after each update.
        pruning_coef: Pruning coefficient for PELT algorithm.
        segment_count: Initial guess for number of segments.
        trim: Trimming proportion for boundary change points.
        momentum_coef: Momentum coefficient for parameter updates.
        multiple_epochs: A function that takes the segment length and
            returns an int for additional epochs.
        epsilon: Epsilon for numerical stability.
        order: Order for AR, VAR, ARIMA, GARCH models.
        p: Number of covariates. If None, inferred from data.
        variance_estimation: Pre-specified variance/covariance matrix.
        cp_only: If True, only change points are returned.
        vanilla_percentage: Interpolation parameter between PELT and SeGD.
        warm_start: If True, use the previous segment's parameters as
            initial values for the new segment.
        **kwargs: Additional model-specific parameters.

    Returns:
        A ``FastCPDResult`` object containing change points, costs,
        residuals, and parameter estimates.
    """
    family = family.lower() if family is not None else 'custom'
    assert family in ('mean', 'variance', 'meanvariance')
    assert cost_adjustment in ('BIC', 'MBIC', 'MDL')

    if variance_estimation is not None:
        variance_estimation = numpy.asarray(variance_estimation)
    elif family == 'mean':
        variance_estimation = fastcpd.variance_estimation.mean(data)

    if family == 'mean':
        p = data.shape[1]
    elif family == 'variance':
        p = data.shape[1] ** 2
    elif family == 'meanvariance':
        p = data.shape[1] + data.shape[1] ** 2

    if pruning_coef is None:
        pruning_coef = 0.0
    if cost_adjustment == "MBIC":
        pruning_coef += p * log(2)
    elif cost_adjustment == "MDL":
        pruning_coef += p * log2(2)

    if isinstance(beta, str):
        if beta == 'BIC':
            beta = (p + 1) * log(data.shape[0]) / 2
        elif beta == 'MBIC':
            beta = (p + 2) * log(data.shape[0]) / 2
        elif beta == 'MDL':
            beta = (p + 2) * log2(data.shape[0]) / 2
        else:
            raise ValueError(f"Unknown beta criterion: {beta}")

    result = fastcpd_impl(
        beta,
        cost_adjustment,
        cp_only,
        data,
        epsilon,
        family,
        line_search,
        [],
        momentum_coef,
        order,
        p,
        0,
        pruning_coef,
        segment_count,
        trim,
        [],
        1.0,
        variance_estimation,
        warm_start,
    )
    return result
