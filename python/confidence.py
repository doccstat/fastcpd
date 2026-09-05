"""
Confidence intervals for fastcpd Python results.
"""

import math
import warnings

import numpy

from ._r_random import RRandom


def _variance_module():
    # Lazy import keeps importing ``fastcpd.confidence`` independent of the
    # package initialisation order.
    from fastcpd import variance_estimation
    return variance_estimation


def confint(
    result,
    data=None,
    parm='cp',
    method=None,
    level=0.95,
    B=999,
    family=None,
    bootstrap='nonparametric',
    window=None,
    min_segment_length=2,
    random_state=None,
    detect_kwargs=None,
    order=None,
):
    """Construct confidence intervals for a ``CpdResult``.

    Args:
        result: A ``fastcpd.segmentation.CpdResult``.
        data: Original data used for fitting. Defaults to ``result.data``.
        parm: ``"cp"`` for change-point locations or ``"theta"`` for segment
            parameters.
        method: For ``parm="cp"``, ``"bootstrap"`` or ``"profile"``. For
            ``parm="theta"``, ``"wald"``. Multivariate linear-model Wald
            intervals are intentionally unsupported until a shared estimator
            contract is defined.
        level: Confidence level.
        B: Number of bootstrap replicates.
        family: Family used to refit bootstrap samples or evaluate profile
            costs. Defaults to ``result.family``. Aliases ``"lm"`` and
            ``"gaussian"`` are equivalent.
        bootstrap: Bootstrap type. Currently only ``"nonparametric"`` is
            implemented.
        window: Optional half-width around each detected change point for
            profile intervals.
        min_segment_length: Minimum observations on each side of a profile
            candidate.
        random_state: Optional NumPy random seed or Generator.
        detect_kwargs: Extra keyword arguments passed to
            ``fastcpd.segmentation.detect`` during bootstrap refits.
        order: Optional model order for legacy result objects that do not
            carry an ``order`` attribute.  For current results the stored
            order is used unless this explicit value is supplied.

    Returns:
        A list of dictionaries. Each dictionary contains the estimate, lower
        and upper interval bounds, and method-specific diagnostics.

    Examples:
        >>> import numpy as np
        >>> from fastcpd import detect_mean
        >>> result = detect_mean(
        ...     np.r_[np.zeros(10), np.full(10, 5.0)],
        ...     beta=2, cost_adjustment='BIC', variance_estimation=np.eye(1)
        ... )
        >>> interval = confint(result, method='profile', level=0.8, window=2)
        >>> interval[0]['estimate']
        10
    """
    if data is None:
        data = getattr(result, 'data', None)
    data = _as_2d_data(data)
    if family is None:
        family = getattr(result, 'family', None)
    if family is None:
        raise ValueError("family must be provided")
    family = str(family).lower()
    result_order = order if order is not None else getattr(result, 'order', None)
    detect_kwargs = _stored_detect_kwargs(
        result, family, result_order, detect_kwargs
    )
    analysis_data = _analysis_data(data, detect_kwargs)

    if not 0 < level < 1:
        raise ValueError("level must be in (0, 1)")
    if parm not in ('cp', 'theta'):
        raise ValueError("parm must be 'cp' or 'theta'")

    if method is None:
        method = 'bootstrap' if parm == 'cp' else 'wald'
    if parm == 'cp' and method == 'bootstrap':
        return _cp_bootstrap(
            result, data, level, B, family, bootstrap, random_state,
            detect_kwargs)
    if parm == 'cp' and method == 'profile':
        return _cp_profile(
            result, analysis_data, level, family, result_order, window,
            min_segment_length)
    if parm == 'theta' and method == 'wald':
        return _theta_wald(result, analysis_data, level, family)

    raise ValueError(f"method {method!r} is not available for parm {parm!r}")


def _as_2d_data(data):
    if data is None:
        raise ValueError("data must be provided")
    data = numpy.asarray(data, dtype=float)
    if data.ndim == 1:
        data = data.reshape(-1, 1)
    if data.ndim != 2 or data.shape[0] == 0 or data.shape[1] == 0:
        raise ValueError("data must be one- or two-dimensional and non-empty")
    if not numpy.all(numpy.isfinite(data)):
        raise ValueError("data must contain only finite numeric values")
    return data


def _stored_detect_kwargs(result, family, order, detect_kwargs):
    """Merge the original fit options with explicit confidence overrides.

    R's confidence implementation refits the complete original call.  The
    Python result stores the equivalent keyword map so bootstrap refits use
    the same penalty, variance estimate, trimming, and model options instead
    of silently falling back to detector defaults.
    """
    options = dict(getattr(result, 'fit_kwargs', {}) or {})
    # Internal metadata keys are useful to wrappers but must never leak to the
    # public ``detect`` function.
    native_family = options.pop('_native_family', None)
    wrapper_family = options.pop('_wrapper_family', None)
    options.pop('data', None)
    options.pop('cp_only', None)
    options.pop('show_progress', None)
    options.update(dict(detect_kwargs or {}))
    # Confidence refits always request change points only and suppress native
    # progress output.  Remove caller-provided values after merging as well;
    # otherwise wrappers that pass ``cp_only=True`` explicitly would receive
    # the same keyword twice from ``_bootstrap_refit``.
    options.pop('cp_only', None)
    options.pop('show_progress', None)
    if native_family is not None:
        options['_native_family'] = native_family
    if wrapper_family is not None:
        options['_wrapper_family'] = wrapper_family
    # ``family``/``order`` supplied to ``confint`` are explicit overrides.
    # This also makes the function usable with pre-metadata result objects;
    # current results simply pass their stored values back through here.
    options['family'] = family
    if order is not None:
        options['order'] = order
    return options


def _analysis_data(data, detect_kwargs):
    """Return the data representation stored by the equivalent R result."""
    if detect_kwargs.get('_wrapper_family') == 'rank':
        # Python retains original observations so bootstrap refits can repeat
        # the public wrapper. Profile and Wald calculations use the centered
        # rank representation consumed by the native mean-family fit.
        from fastcpd.segmentation import _rank_transform
        return _rank_transform(data)
    return data


def _normalize_family(family):
    if family is None:
        raise ValueError("family must be provided")
    family = str(family).lower()
    if family == 'gaussian':
        return 'lm'
    if family == 'mgaussian':
        return 'var'
    if family == 'rank':
        return 'mean'
    return family


def _rng(random_state):
    if isinstance(random_state, (
        RRandom, numpy.random.Generator, numpy.random.RandomState
    )):
        return random_state
    if isinstance(random_state, (int, numpy.integer)):
        return RRandom(random_state)
    return numpy.random.default_rng(random_state)


def _cp_bootstrap(
    result, data, level, B, family, bootstrap, random_state, detect_kwargs
):
    if bootstrap != 'nonparametric':
        raise NotImplementedError(
            "Only bootstrap='nonparametric' is currently implemented")
    if family is None:
        raise ValueError("family must be provided")
    family = str(family).lower()
    data = _as_2d_data(data)
    B = int(B)
    if B <= 0:
        raise ValueError("B must be a positive integer")
    rng = _rng(random_state)
    reference_cp = sorted(int(cp) for cp in result.cp_set)
    if not reference_cp:
        return []

    matched = numpy.full((B, len(reference_cp)), numpy.nan)
    from fastcpd import segmentation

    failed_refits = 0
    for b in range(B):
        boot_data = _segment_bootstrap_data(data, reference_cp, rng)
        try:
            boot_result = _bootstrap_refit(
                segmentation, boot_data, family, detect_kwargs, rng
            )
            boot_cp = boot_result.cp_set
        except Exception:
            failed_refits += 1
            boot_cp = []
        matched[b, :] = _match_cp_set(
            reference_cp, sorted(int(cp) for cp in boot_cp), data.shape[0])

    alpha = 1 - level
    rows = []
    for i, estimate in enumerate(reference_cp):
        estimates = matched[:, i]
        detected = estimates[~numpy.isnan(estimates)]
        if detected.size:
            lower, upper = _quantile_type1(
                detected, (alpha / 2, 1 - alpha / 2)
            )
        else:
            lower, upper = math.nan, math.nan
        rows.append({
            'parm': 'cp',
            'index': i + 1,
            'estimate': estimate,
            'lower': float(lower),
            'upper': float(upper),
            'detection_rate': float(detected.size / B),
            'level': level,
            'method': 'bootstrap',
            'bootstrap': bootstrap,
        })
    if failed_refits:
        warnings.warn(
            f"{failed_refits} bootstrap refit(s) failed; their intervals "
            "use NaN.",
            RuntimeWarning,
            stacklevel=3,
        )
    return rows


def _bootstrap_refit(segmentation, data, family, detect_kwargs, rng):
    """Refit a bootstrap sample using the stored wrapper contract."""
    options = dict(detect_kwargs)
    wrapper_family = options.pop('_wrapper_family', None)
    native_family = options.pop('_native_family', None)
    # These are controlled by the bootstrap driver.  In particular, popping
    # them here protects compatibility with old ``fit_kwargs`` maps that may
    # still contain ``cp_only``/``show_progress`` from a prior release.
    options.pop('cp_only', None)
    options.pop('show_progress', None)
    fit_family = wrapper_family or options.get('family') or family
    # The wrapper marker is authoritative when a result was produced through
    # rank/KCP; those transforms must be repeated on each bootstrap sample.
    if fit_family == 'rank':
        options.pop('family', None)
        return segmentation.detect_rank(
            data, cp_only=True, **options
        )
    if fit_family in ('kernel', 'kcp'):
        options.pop('family', None)
        # R uses one advancing stream for segment resampling and every KCP
        # refit. Scalar Python seeds follow that stream exactly; native NumPy
        # generators retain the historical per-refit seed behavior.
        options['random_state'] = (
            rng if isinstance(rng, RRandom) else _draw_seed(rng)
        )
        return segmentation.detect_kcp(
            data, cp_only=True, **options
        )
    options['family'] = native_family or options.get('family') or family
    options['cp_only'] = True
    options['show_progress'] = False
    return segmentation.detect(data=data, **options)


def _draw_seed(rng):
    """Draw a uint32-compatible seed from either NumPy RNG API."""
    if hasattr(rng, 'integers'):
        return int(rng.integers(0, 2**32 - 1))
    # Windows gives RandomState's default C ``long`` dtype only 32 signed
    # bits, making the uint32 upper bound invalid. Request uint32 explicitly;
    # this preserves the same random draw on every platform.
    return int(rng.randint(0, 2**32 - 1, dtype=numpy.uint32))


def _segment_bootstrap_data(data, cp_set, rng):
    bounds = [0] + list(cp_set) + [data.shape[0]]
    boot_data = numpy.array(data, copy=True)
    for start, end in zip(bounds[:-1], bounds[1:]):
        rows = numpy.arange(start, end)
        if rows.size:
            boot_data[rows, :] = data[rng.choice(rows, rows.size, replace=True), :]
    return boot_data


def _match_cp_set(reference_cp, bootstrap_cp, n):
    if not bootstrap_cp:
        return [math.nan] * len(reference_cp)
    matched = []
    for i, cp in enumerate(reference_cp):
        left = 0 if i == 0 else math.floor((reference_cp[i - 1] + cp) / 2)
        right = n if i == len(reference_cp) - 1 else math.ceil(
            (cp + reference_cp[i + 1]) / 2)
        candidates = [x for x in bootstrap_cp if left < x <= right]
        if candidates:
            matched.append(min(candidates, key=lambda x: abs(x - cp)))
        else:
            matched.append(math.nan)
    return matched


def _quantile_type1(values, probabilities):
    """R's ``stats::quantile(..., type = 1)`` for finite 1-D values."""
    values = numpy.sort(numpy.asarray(values, dtype=float))
    probabilities = tuple(probabilities)
    if values.size == 0:
        return (math.nan,) * len(probabilities)
    out = []
    n = values.size
    for probability in probabilities:
        probability = float(probability)
        if probability <= 0:
            index = 0
        elif probability >= 1:
            index = n - 1
        else:
            # Inverse empirical CDF: min{i/n >= p}, with one-based i.
            index = int(math.ceil(probability * n)) - 1
            index = max(0, min(n - 1, index))
        out.append(float(values[index]))
    return tuple(out)


def _cp_profile(
    result, data, level, family, order, window, min_segment_length
):
    family = _normalize_family(family)
    data = _as_2d_data(data)
    fit_options = dict(getattr(result, 'fit_kwargs', {}) or {})
    # ``confint()`` has already converted rank-wrapper input to centered ranks
    # through _analysis_data(). Keep the mean-family dispatch marker here
    # without applying that transform a second time.
    if str(fit_options.get('family', '')).lower() == 'rank':
        family = 'mean'
    variance_estimation = fit_options.get('variance_estimation')
    cost = _profile_cost_function(
        data, family, order, variance_estimation=variance_estimation
    )
    cp_set = sorted(int(cp) for cp in result.cp_set)
    if not cp_set:
        return []

    try:
        min_segment_length_float = float(min_segment_length)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("min_segment_length must be a positive integer") from error
    min_segment_length = int(min_segment_length_float)
    if (min_segment_length_float != min_segment_length or
            min_segment_length <= 0):
        raise ValueError("min_segment_length must be positive")
    if window is not None:
        try:
            window_float = float(window)
        except (TypeError, ValueError, OverflowError) as error:
            raise ValueError("window must be a non-negative integer") from error
        if window_float < 0 or window_float != int(window_float):
            raise ValueError("window must be a non-negative integer")
        window = int(window_float)
    cutoff = _chisq1(level) / 2
    bounds = [0] + cp_set + [data.shape[0]]
    rows = []

    for i, cp in enumerate(cp_set):
        left = bounds[i]
        right = bounds[i + 2]
        tau_min = left + min_segment_length
        tau_max = right - min_segment_length
        if window is not None:
            tau_min = max(tau_min, cp - window)
            tau_max = min(tau_max, cp + window)
        candidates = list(range(tau_min, tau_max + 1))
        costs = [
            cost(left, tau) + cost(tau, right)
            for tau in candidates
        ]
        finite = [(tau, value) for tau, value in zip(candidates, costs)
                  if math.isfinite(value)]
        if finite:
            profile_min = min(value for _, value in finite)
            support = [
                tau for tau, value in finite if value - profile_min <= cutoff
            ]
            lower = min(min(support), cp)
            upper = max(max(support), cp)
        else:
            profile_min, lower, upper = math.nan, math.nan, math.nan
        rows.append({
            'parm': 'cp',
            'index': i + 1,
            'estimate': cp,
            'lower': lower,
            'upper': upper,
            'profile_min': profile_min,
            'cutoff': cutoff,
            'level': level,
            'method': 'profile',
        })
    return rows


def _profile_cost_function(data, family, order=None, variance_estimation=None):
    if family == 'mean':
        if variance_estimation is None:
            try:
                variance_estimation = (
                    _variance_module().estimate_variance_mean(data)
                )
            except Exception:
                variance_estimation = numpy.eye(data.shape[1])
        sigma_inv = _precision(variance_estimation, data.shape[1])

        def cost(start, end):
            segment = data[start:end, :]
            centered = segment - segment.mean(axis=0)
            return float(numpy.sum((centered @ sigma_inv) * centered) / 2)
        return cost

    if family == 'variance':
        centered_data = data - data.mean(axis=0)

        def cost(start, end):
            segment = centered_data[start:end, :]
            covariance = segment.T @ segment / segment.shape[0]
            return float(segment.shape[0] * _logdet(covariance) / 2)
        return cost

    if family == 'meanvariance':
        def cost(start, end):
            segment = data[start:end, :]
            centered = segment - segment.mean(axis=0)
            covariance = centered.T @ centered / segment.shape[0]
            return float(segment.shape[0] * _logdet(covariance) / 2)
        return cost

    if family == 'exponential':
        def cost(start, end):
            segment = data[start:end, 0]
            if numpy.any(segment <= 0):
                return math.inf
            return float(segment.size * (numpy.log(segment.mean()) + 1))
        return cost

    if family in ('lm', 'gaussian'):
        sigma2 = _linear_regression_variance(data, variance_estimation)

        def cost(start, end):
            segment = data[start:end, :]
            x = segment[:, 1:]
            y = segment[:, 0]
            if x.shape[0] <= x.shape[1]:
                return math.inf
            try:
                beta, _, _, _ = numpy.linalg.lstsq(x, y, rcond=None)
            except numpy.linalg.LinAlgError:
                return math.inf
            residual = y - x @ beta
            return float(residual @ residual / (2 * sigma2))
        return cost

    if family in ('binomial', 'poisson'):
        def cost(start, end):
            segment = data[start:end, :]
            x = segment[:, 1:]
            if x.shape[0] <= x.shape[1]:
                return math.inf
            return _native_single_segment_cost(
                segment, family=family, order=(0, 0, 0)
            )
        return cost

    if family == 'quantile':
        if order is None or len(order) < 1:
            raise ValueError("quantile profile intervals require result.order")
        tau = float(order[0])

        def cost(start, end):
            segment = data[start:end, :]
            x = segment[:, 1:]
            if x.shape[0] <= x.shape[1]:
                return math.inf
            return _native_single_segment_cost(
                segment, family='quantile', order=(tau,)
            )
        return cost

    if family == 'arima':
        if order is None or len(order) != 3:
            raise ValueError("ARIMA profile intervals require result.order")
        arima_order = tuple(int(value) for value in order)

        def cost(start, end):
            return _native_single_segment_cost(
                data[start:end, 0], family='arima', order=arima_order
            )
        return cost

    raise NotImplementedError(
        f"Profile intervals are not implemented for family {family!r}")


def _precision(variance_estimation, dimension):
    """Return an R-compatible precision matrix, with identity fallback."""
    if variance_estimation is None:
        return numpy.eye(dimension)
    sigma = numpy.asarray(variance_estimation, dtype=float)
    if sigma.ndim == 0:
        sigma = numpy.array([[float(sigma)]])
    if sigma.ndim != 2 or sigma.shape[0] != sigma.shape[1]:
        return numpy.eye(dimension)
    if sigma.shape != (dimension, dimension) or not numpy.all(
            numpy.isfinite(sigma)):
        return numpy.eye(dimension)
    try:
        return numpy.linalg.inv(sigma)
    except numpy.linalg.LinAlgError:
        return numpy.eye(dimension)


def _linear_regression_variance(data, supplied=None):
    if supplied is None:
        try:
            value = _variance_module().estimate_variance_linear_regression(data)
        except Exception:
            value = 1.0
    else:
        values = numpy.asarray(supplied, dtype=float).ravel()
        value = values[0] if values.size else 1.0
    try:
        value = float(value)
    except (TypeError, ValueError, OverflowError):
        value = 1.0
    return value if math.isfinite(value) and value > 0 else 1.0


def _theta_wald(result, data, level, family):
    family = _normalize_family(family)
    data = _as_2d_data(data)
    if bool(getattr(result, 'cp_only', False)):
        raise ValueError(
            "Wald intervals require a result fitted with cp_only=False"
        )
    theta = numpy.asarray(result.thetas, dtype=float)
    if theta.ndim == 1:
        theta = theta.reshape(1, -1)
    if theta.size == 0:
        raise ValueError("Wald intervals require a result with thetas")

    fit_options = dict(getattr(result, 'fit_kwargs', {}) or {})
    p_response = fit_options.get('p_response', 0)
    try:
        p_response = int(p_response)
    except (TypeError, ValueError, OverflowError):
        p_response = 0
    if family in ('lm', 'gaussian') and p_response > 1:
        raise NotImplementedError(
            "Wald intervals are not implemented for multivariate LM results"
        )
    se = _theta_se_function(
        data, family, p_response=p_response,
        parameter_count=int(theta.shape[0]),
    )
    z_value = _normal_quantile(1 - (1 - level) / 2)
    bounds = [0] + sorted(int(cp) for cp in result.cp_set) + [data.shape[0]]
    if theta.shape[1] != len(bounds) - 1:
        raise ValueError(
            "result.thetas does not have one column per detected segment"
        )
    rows = []
    for segment_index, (start, end) in enumerate(zip(bounds[:-1], bounds[1:])):
        segment_theta = theta[:, segment_index]
        segment_se = se(start, end, segment_theta)
        for param_index, estimate in enumerate(theta[:, segment_index]):
            se_value = segment_se[param_index]
            rows.append({
                'parm': 'theta',
                'segment': segment_index + 1,
                'parameter': param_index + 1,
                'estimate': float(estimate),
                'lower': float(estimate - z_value * se_value),
                'upper': float(estimate + z_value * se_value),
                'se': float(se_value),
                'level': level,
                'method': 'wald',
            })
    return rows


def _theta_se_function(data, family, *, p_response=0, parameter_count=None):
    if family == 'mean':
        def se(start, end, estimate=None):
            segment = data[start:end, :]
            if segment.shape[0] <= 1:
                return numpy.full(segment.shape[1], numpy.nan)
            covariance = numpy.atleast_2d(numpy.cov(segment, rowvar=False))
            return numpy.sqrt(numpy.diag(covariance) /
                              segment.shape[0])
        return se

    if family == 'exponential':
        def se(start, end, estimate=None):
            segment = data[start:end, 0]
            rate = 1 / segment.mean()
            return numpy.array([rate / numpy.sqrt(segment.size)])
        return se

    if family in ('lm', 'gaussian'):
        # ``lm(..., p_response=q)`` is dispatched to the native multivariate
        # Gaussian family.  Its flattened coefficient vector follows R/
        # Armadillo's column-major response-block layout (q blocks of k
        # predictor coefficients).  Infer q when older manually constructed
        # results do not carry fit metadata, then use the usual
        # Sigma ⊗ (X'X)^-1 covariance diagonal.
        response_count = _infer_response_count(
            data.shape[1], p_response, parameter_count
        )
        if response_count > 1:
            def se(start, end, estimate=None):
                segment = data[start:end, :]
                q = response_count
                x = segment[:, q:]
                y = segment[:, :q]
                predictor_count = x.shape[1]
                if (predictor_count == 0 or
                        segment.shape[0] <= predictor_count):
                    return numpy.full(q * predictor_count, numpy.nan)
                try:
                    xtx_inv = numpy.linalg.inv(x.T @ x)
                    beta, _, _, _ = numpy.linalg.lstsq(x, y, rcond=None)
                except numpy.linalg.LinAlgError:
                    return numpy.full(q * predictor_count, numpy.nan)
                residual = y - x @ beta
                degrees = max(segment.shape[0] - predictor_count, 1)
                sigma = (residual.T @ residual) / degrees
                diagonal = (
                    numpy.diag(sigma)[:, None] *
                    numpy.diag(xtx_inv)[None, :]
                )
                diagonal = diagonal.reshape(-1, order='C')
                if (diagonal.size != q * predictor_count or
                        numpy.any(diagonal < 0) or
                        not numpy.all(numpy.isfinite(diagonal))):
                    return numpy.full(q * predictor_count, numpy.nan)
                return numpy.sqrt(diagonal)
            return se

        def se(start, end, estimate=None):
            segment = data[start:end, :]
            x = segment[:, 1:]
            y = segment[:, 0]
            try:
                beta, _, _, _ = numpy.linalg.lstsq(x, y, rcond=None)
                residual = y - x @ beta
                sigma2 = residual @ residual / max(x.shape[0] - x.shape[1], 1)
                xtx_inv = numpy.linalg.inv(x.T @ x)
            except numpy.linalg.LinAlgError:
                return numpy.full(x.shape[1], numpy.nan)
            return numpy.sqrt(numpy.diag(xtx_inv) * sigma2)
        return se

    if family in ('binomial', 'poisson'):
        def se(start, end, estimate=None):
            segment = data[start:end, :]
            x = segment[:, 1:]
            try:
                # R's Wald method refits a fresh IRLS model for each segment;
                # using the detection theta can give stale/zero weights after
                # a SeGD fit.  The native one-segment path shares the same
                # likelihood and provides a stable fresh coefficient vector.
                fitted = _native_single_segment_fit(
                    segment, family=family, order=(0, 0, 0)
                )
                if fitted is None:
                    return numpy.full(x.shape[1], numpy.nan)
                fitted_theta = numpy.asarray(fitted.thetas, dtype=float)
                if fitted_theta.ndim == 2:
                    fitted_theta = fitted_theta[:, 0]
                fitted_theta = fitted_theta.reshape(-1)
                if fitted_theta.size != x.shape[1]:
                    return numpy.full(x.shape[1], numpy.nan)
                # Keep the vectors one-dimensional.  A 1-column matrix from
                # the native result otherwise broadcasts against ``weights``
                # and produces an n×n information matrix.
                eta = numpy.asarray(x @ fitted_theta, dtype=float).reshape(-1)
                if family == 'binomial':
                    mu = _logistic(eta)
                    weights = mu * (1 - mu)
                else:
                    with numpy.errstate(over='raise', invalid='raise'):
                        weights = numpy.exp(eta)
                weights = numpy.asarray(weights, dtype=float).reshape(-1)
                if weights.size != x.shape[0] or not numpy.all(
                        numpy.isfinite(weights)):
                    return numpy.full(x.shape[1], numpy.nan)
                separation_threshold = math.sqrt(numpy.finfo(float).eps)
                if numpy.any(weights <= separation_threshold):
                    return numpy.full(x.shape[1], numpy.nan)
                information = x.T @ (x * weights[:, None])
                condition = numpy.linalg.cond(information, p=1)
                if (not math.isfinite(condition) or
                        condition >= 1 / separation_threshold):
                    return numpy.full(x.shape[1], numpy.nan)
                information_inv = numpy.linalg.inv(information)
            except (FloatingPointError, ValueError, numpy.linalg.LinAlgError,
                    RuntimeError):
                return numpy.full(x.shape[1], numpy.nan)
            diagonal = numpy.diag(information_inv)
            if numpy.any(diagonal < 0) or not numpy.all(
                    numpy.isfinite(diagonal)):
                return numpy.full(x.shape[1], numpy.nan)
            return numpy.sqrt(diagonal)
        return se

    raise NotImplementedError(
        f"Wald intervals are not implemented for family {family!r}")


def _infer_response_count(column_count, p_response=0, parameter_count=None):
    """Infer the response block size for a flattened LM result.

    A multivariate regression data matrix has ``q`` response columns and
    ``k`` predictors, while native ``thetas`` has ``q*k`` entries.  The
    explicit ``p_response`` metadata is authoritative; for legacy objects we
    solve the small integer factorisation and choose the smaller response
    block when the factorisation is symmetric/ambiguous.
    """
    try:
        q = int(p_response)
    except (TypeError, ValueError, OverflowError):
        q = 0
    if q > 0 and q < column_count:
        if parameter_count is None or q * (column_count - q) == parameter_count:
            return q
    if parameter_count is None:
        return 1
    candidates = [
        value for value in range(1, column_count)
        if value * (column_count - value) == parameter_count
    ]
    if not candidates:
        return 1
    return min(candidates)


def _logistic(eta):
    """Evaluate the inverse-logit without overflowing either tail."""
    eta = numpy.asarray(eta, dtype=float)
    result = numpy.empty_like(eta)
    positive = eta >= 0
    result[positive] = 1 / (1 + numpy.exp(-eta[positive]))
    exp_eta = numpy.exp(eta[~positive])
    result[~positive] = exp_eta / (1 + exp_eta)
    return result


def _native_single_segment_cost(data, family, order):
    """Evaluate one segment with the shared native R/Python likelihood."""
    fit = _native_single_segment_fit(data, family, order)
    if fit is None or fit.cp_set.size or fit.cost_values.size != 1:
        return math.inf
    return float(fit.cost_values[0])


def _native_single_segment_fit(data, family, order):
    """Fit one segment through the native detector, or return ``None``."""
    from fastcpd import segmentation
    try:
        fit = segmentation.detect(
            data=numpy.asarray(data, dtype=float),
            family=family,
            order=order,
            beta=1e100,
            cost_adjustment='BIC',
            segment_count=1,
            trim=0.0,
            vanilla_percentage=1.0,
            cp_only=False,
            show_progress=False,
        )
    except (RuntimeError, ValueError, numpy.linalg.LinAlgError,
            FloatingPointError):
        return None
    return fit


def _logdet(matrix):
    sign, value = numpy.linalg.slogdet(matrix)
    if sign <= 0 or not math.isfinite(value):
        values = numpy.linalg.eigvalsh(matrix)
        value = numpy.sum(numpy.log(numpy.maximum(values, numpy.finfo(float).eps)))
    return value


def _normal_quantile(p):
    return _inverse_standard_normal(p)


def _chisq1(p):
    z = _inverse_standard_normal((1 + p) / 2)
    return z * z


def _inverse_standard_normal(p):
    # Peter J. Acklam's rational approximation followed by one Halley
    # correction against the standard-normal CDF.  The correction brings the
    # common confidence levels to the same double-precision values returned
    # by R's qnorm(), which also keeps profile cutoffs numerically aligned.
    if not 0 < p < 1:
        raise ValueError("p must be in (0, 1)")
    a = [
        -3.969683028665376e+01, 2.209460984245205e+02,
        -2.759285104469687e+02, 1.383577518672690e+02,
        -3.066479806614716e+01, 2.506628277459239e+00,
    ]
    b = [
        -5.447609879822406e+01, 1.615858368580409e+02,
        -1.556989798598866e+02, 6.680131188771972e+01,
        -1.328068155288572e+01,
    ]
    c = [
        -7.784894002430293e-03, -3.223964580411365e-01,
        -2.400758277161838e+00, -2.549732539343734e+00,
        4.374664141464968e+00, 2.938163982698783e+00,
    ]
    d = [
        7.784695709041462e-03, 3.224671290700398e-01,
        2.445134137142996e+00, 3.754408661907416e+00,
    ]
    plow = 0.02425
    phigh = 1 - plow
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        estimate = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q +
                     c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) *
                                            q + d[3]) * q + 1)
    elif p <= phigh:
        q = p - 0.5
        r = q * q
        estimate = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r +
                     a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r +
                                                 b[2]) * r + b[3]) * r +
                                               b[4]) * r + 1)
    else:
        q = math.sqrt(-2 * math.log(1 - p))
        estimate = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q +
                      c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) *
                                            q + d[3]) * q + 1)

    cdf = 0.5 * math.erfc(-estimate / math.sqrt(2.0))
    density = math.exp(-estimate * estimate / 2.0) / math.sqrt(2.0 * math.pi)
    if density > 0:
        correction = (cdf - p) / density
        estimate -= correction / (1.0 + estimate * correction / 2.0)
    return estimate
