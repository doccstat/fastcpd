import numpy
from numpy.typing import NDArray


def fastcpd_impl(
    beta: float,
    cost_adjustment: str,
    cp_only: bool,
    data: NDArray[numpy.float64],
    epsilon: float,
    family: str,
    line_search: NDArray[numpy.float64],
    lower: NDArray[numpy.float64],
    momentum_coef: float,
    order: NDArray[numpy.int64],
    p: int,
    p_response: int,
    pruning_coef: float,
    segment_count: int,
    trim: float,
    upper: NDArray[numpy.float64],
    vanilla_percentage: float,
    variance_estimate: NDArray[numpy.float64],
    warm_start: bool,
) -> NDArray[numpy.float64]:
    ...
