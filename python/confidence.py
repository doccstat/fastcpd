"""
Confidence intervals for fastcpd Python results.
"""

import math

import numpy


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
):
    """Construct confidence intervals for a ``CpdResult``.

    Args:
        result: A ``fastcpd.segmentation.CpdResult``.
        data: Original data used for fitting. Required because Python result
            objects intentionally keep only compact fit output.
        parm: ``"cp"`` for change-point locations or ``"theta"`` for segment
            parameters.
        method: For ``parm="cp"``, ``"bootstrap"`` or ``"profile"``. For
            ``parm="theta"``, ``"wald"``.
        level: Confidence level.
        B: Number of bootstrap replicates.
        family: Family used to refit bootstrap samples or evaluate profile
            costs. Aliases ``"lm"`` and ``"gaussian"`` are equivalent.
        bootstrap: Bootstrap type. Currently only ``"nonparametric"`` is
            implemented.
        window: Optional half-width around each detected change point for
            profile intervals.
        min_segment_length: Minimum observations on each side of a profile
            candidate.
        random_state: Optional NumPy random seed or Generator.
        detect_kwargs: Extra keyword arguments passed to
            ``fastcpd.segmentation.detect`` during bootstrap refits.

    Returns:
        A list of dictionaries. Each dictionary contains the estimate, lower
        and upper interval bounds, and method-specific diagnostics.
    """
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
            result, data, level, family, window, min_segment_length)
    if parm == 'theta' and method == 'wald':
        return _theta_wald(result, data, level, family)

    raise ValueError(f"method {method!r} is not available for parm {parm!r}")


def _as_2d_data(data):
    if data is None:
        raise ValueError("data must be provided")
    data = numpy.asarray(data, dtype=float)
    if data.ndim == 1:
        data = data.reshape(-1, 1)
    return data


def _normalize_family(family):
    if family is None:
        raise ValueError("family must be provided")
    family = family.lower()
    if family == 'lm':
        return 'gaussian'
    if family == 'var':
        return 'mgaussian'
    return family


def _rng(random_state):
    if isinstance(random_state, numpy.random.Generator):
        return random_state
    return numpy.random.default_rng(random_state)


def _cp_bootstrap(
    result, data, level, B, family, bootstrap, random_state, detect_kwargs
):
    if bootstrap != 'nonparametric':
        raise NotImplementedError(
            "Only bootstrap='nonparametric' is currently implemented")
    family = _normalize_family(family)
    data = _as_2d_data(data)
    B = int(B)
    if B <= 0:
        raise ValueError("B must be a positive integer")
    detect_kwargs = {} if detect_kwargs is None else dict(detect_kwargs)
    rng = _rng(random_state)
    reference_cp = sorted(int(cp) for cp in result.cp_set)
    if not reference_cp:
        return []

    matched = numpy.full((B, len(reference_cp)), numpy.nan)
    from fastcpd import segmentation

    for b in range(B):
        boot_data = _segment_bootstrap_data(data, reference_cp, rng)
        try:
            boot_cp = segmentation.detect(
                data=boot_data,
                family=family,
                cp_only=True,
                **detect_kwargs,
            )
        except Exception:
            boot_cp = []
        matched[b, :] = _match_cp_set(
            reference_cp, sorted(int(cp) for cp in boot_cp), data.shape[0])

    alpha = 1 - level
    rows = []
    for i, estimate in enumerate(reference_cp):
        estimates = matched[:, i]
        detected = estimates[~numpy.isnan(estimates)]
        if detected.size:
            lower, upper = numpy.quantile(
                detected, [alpha / 2, 1 - alpha / 2])
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
    return rows


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


def _cp_profile(result, data, level, family, window, min_segment_length):
    family = _normalize_family(family)
    data = _as_2d_data(data)
    cost = _profile_cost_function(data, family)
    cp_set = sorted(int(cp) for cp in result.cp_set)
    if not cp_set:
        return []

    min_segment_length = int(min_segment_length)
    if min_segment_length <= 0:
        raise ValueError("min_segment_length must be positive")
    cutoff = _chisq1(level) / 2
    bounds = [0] + cp_set + [data.shape[0]]
    rows = []

    for i, cp in enumerate(cp_set):
        left = bounds[i]
        right = bounds[i + 2]
        tau_min = left + min_segment_length
        tau_max = right - min_segment_length
        if window is not None:
            tau_min = max(tau_min, cp - int(window))
            tau_max = min(tau_max, cp + int(window))
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
            lower, upper = min(support), max(support)
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


def _profile_cost_function(data, family):
    if family == 'mean':
        def cost(start, end):
            segment = data[start:end, :]
            centered = segment - segment.mean(axis=0)
            return float(numpy.sum(centered * centered) / 2)
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

    if family == 'gaussian':
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
            return float(residual @ residual / 2)
        return cost

    raise NotImplementedError(
        f"Profile intervals are not implemented for family {family!r}")


def _theta_wald(result, data, level, family):
    family = _normalize_family(family)
    data = _as_2d_data(data)
    theta = numpy.asarray(result.thetas, dtype=float)
    if theta.ndim == 1:
        theta = theta.reshape(1, -1)
    if theta.size == 0:
        raise ValueError("Wald intervals require a result with thetas")

    se = _theta_se_function(data, family)
    z_value = _normal_quantile(1 - (1 - level) / 2)
    bounds = [0] + sorted(int(cp) for cp in result.cp_set) + [data.shape[0]]
    rows = []
    for segment_index, (start, end) in enumerate(zip(bounds[:-1], bounds[1:])):
        segment_se = se(start, end)
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


def _theta_se_function(data, family):
    if family == 'mean':
        def se(start, end):
            segment = data[start:end, :]
            if segment.shape[0] <= 1:
                return numpy.full(segment.shape[1], numpy.nan)
            covariance = numpy.atleast_2d(numpy.cov(segment, rowvar=False))
            return numpy.sqrt(numpy.diag(covariance) /
                              segment.shape[0])
        return se

    if family == 'exponential':
        def se(start, end):
            segment = data[start:end, 0]
            rate = 1 / segment.mean()
            return numpy.array([rate / numpy.sqrt(segment.size)])
        return se

    if family == 'gaussian':
        def se(start, end):
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

    raise NotImplementedError(
        f"Wald intervals are not implemented for family {family!r}")


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
    # Peter J. Acklam's rational approximation, sufficient for CI endpoints.
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
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q +
                 c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) *
                                        q + d[3]) * q + 1)
    if p <= phigh:
        q = p - 0.5
        r = q * q
        return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r +
                 a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r +
                                             b[2]) * r + b[3]) * r +
                                           b[4]) * r + 1)
    q = math.sqrt(-2 * math.log(1 - p))
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q +
              c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) *
                                    q + d[3]) * q + 1)
