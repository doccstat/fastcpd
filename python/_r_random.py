"""R-compatible random primitives for deterministic cross-language results.

Only scalar integer seeds use this implementation. NumPy ``Generator`` and
``RandomState`` inputs retain their native streams. The generator installs
R's historical Mersenne-Twister seed state into NumPy's compiled MT19937 bit
generator, so bulk uniform draws stay vectorized and stochastic preprocessing
does not enter the native fastcpd detector loop.
"""

from __future__ import annotations

from math import ceil, log2, prod
from statistics import NormalDist
from typing import TypeAlias

import numpy


Seed: TypeAlias = int | numpy.integer
Shape: TypeAlias = int | tuple[int, ...]

_UINT32_MASK = (1 << 32) - 1
_UINT32_SCALE = 2.0**-32
_R_FIXUP = 0.5 / ((1 << 32) - 1)
_R_NORMAL_PRECISION = 1 << 27
_STANDARD_NORMAL = NormalDist()


def is_r_seed(value: object) -> bool:
    """Return whether *value* is a scalar integer seed accepted by R."""
    if isinstance(value, (bool, numpy.bool_)):
        return False
    if not isinstance(value, (int, numpy.integer)):
        return False
    seed = int(value)
    return -(1 << 31) < seed < (1 << 31)


class RRandom:
    """Subset of R's default RNG used by KCP and bootstrap resampling."""

    def __init__(self, seed: Seed) -> None:
        if not is_r_seed(seed):
            raise ValueError(
                "R-compatible seeds must be integers between "
                "-2147483647 and 2147483647"
            )
        self._bit_generator = numpy.random.MT19937()
        state = self._bit_generator.state
        state["state"]["key"] = _r_mt_key(int(seed))
        state["state"]["pos"] = 624
        self._bit_generator.state = state

    def random(self, size: Shape | None = None) -> float | numpy.ndarray:
        """Draw R ``runif(0, 1)`` values in generation order."""
        shape, count = _shape_and_count(size)
        values = self._raw(count).astype(numpy.float64) * _UINT32_SCALE
        zero = values == 0.0
        if numpy.any(zero):
            values[zero] = _R_FIXUP
        if shape is None:
            return float(values[0])
        return values.reshape(shape)

    def uniform(
        self,
        low: float = 0.0,
        high: float = 1.0,
        size: Shape | None = None,
    ) -> float | numpy.ndarray:
        """Draw values with R's ``runif`` scaling."""
        values = self.random(size)
        return low + (high - low) * values

    def normal(
        self,
        loc: float = 0.0,
        scale: float = 1.0,
        size: Shape | None = None,
    ) -> float | numpy.ndarray:
        """Draw values with R's default inversion-normal convention."""
        if scale < 0 or not numpy.isfinite(scale):
            raise ValueError("normal scale must be finite and non-negative")
        shape, count = _shape_and_count(size)
        uniforms = numpy.asarray(self.random(2 * count)).reshape(count, 2)
        probabilities = (
            numpy.floor(_R_NORMAL_PRECISION * uniforms[:, 0])
            + uniforms[:, 1]
        ) / _R_NORMAL_PRECISION
        values = numpy.fromiter(
            (_STANDARD_NORMAL.inv_cdf(float(p)) for p in probabilities),
            dtype=numpy.float64,
            count=count,
        )
        values = loc + scale * values
        if shape is None:
            return float(values[0])
        # R matrices and arrays are filled in column-major generation order.
        return values.reshape(shape, order="F")

    def choice(
        self,
        population: int | numpy.ndarray,
        size: int,
        replace: bool = True,
    ) -> numpy.ndarray:
        """Sample uniformly with R 3.6+'s rejection-sampling convention."""
        size = int(size)
        if size < 0:
            raise ValueError("sample size must be non-negative")
        if isinstance(population, (int, numpy.integer)):
            n = int(population)
            values: numpy.ndarray | None = None
        else:
            values = numpy.asarray(population)
            if values.ndim != 1:
                raise ValueError("sample population must be one-dimensional")
            n = int(values.size)
        if n <= 0 and size:
            raise ValueError("cannot sample from an empty population")
        if not replace and size > n:
            raise ValueError("sample size exceeds the population")

        if replace:
            indices = self._unif_indices(n, size)
        else:
            available = numpy.arange(n, dtype=numpy.int64)
            indices = numpy.empty(size, dtype=numpy.int64)
            remaining = n
            for offset in range(size):
                selected = int(self._unif_indices(remaining, 1)[0])
                indices[offset] = available[selected]
                remaining -= 1
                available[selected] = available[remaining]
        return indices if values is None else values[indices]

    def _raw(self, count: int) -> numpy.ndarray:
        return numpy.asarray(
            self._bit_generator.random_raw(count), dtype=numpy.uint64
        )

    def _unif_indices(self, population_size: int, count: int) -> numpy.ndarray:
        if count == 0:
            return numpy.empty(0, dtype=numpy.int64)
        bits = int(ceil(log2(population_size))) if population_size > 1 else 0
        chunks = bits // 16 + 1
        mask = (1 << bits) - 1
        accepted: list[numpy.ndarray] = []
        missing = count
        while missing:
            raw = self._raw(missing * chunks).reshape(missing, chunks)
            words = raw >> 16
            candidates = numpy.zeros(missing, dtype=numpy.uint64)
            for column in range(chunks):
                candidates = (candidates << 16) + words[:, column]
            candidates &= numpy.uint64(mask)
            valid = candidates[candidates < population_size]
            if valid.size:
                accepted.append(valid.astype(numpy.int64, copy=False))
                missing -= int(valid.size)
        return numpy.concatenate(accepted)


def _shape_and_count(size: Shape | None) -> tuple[tuple[int, ...] | None, int]:
    if size is None:
        return None, 1
    if isinstance(size, (int, numpy.integer)):
        shape = (int(size),)
    else:
        shape = tuple(int(dimension) for dimension in size)
    if any(dimension < 0 for dimension in shape):
        raise ValueError("negative dimensions are not allowed")
    return shape, prod(shape)


def _r_mt_key(seed: int) -> numpy.ndarray:
    """Construct the 624-word MT19937 state produced by R ``set.seed``."""
    state = seed & _UINT32_MASK
    for _ in range(50):
        state = (69069 * state + 1) & _UINT32_MASK
    words = numpy.empty(625, dtype=numpy.uint32)
    for index in range(625):
        state = (69069 * state + 1) & _UINT32_MASK
        words[index] = state
    return words[1:]
