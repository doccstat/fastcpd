## Updates since last CRAN release (0.10.3)

### fastcpd 0.11.0

*   Add penalty selection criteria using

    1. BIC: `(p + 1) * log(nrow(data)) / 2`
    1. Modified BIC: `(p + 2) * log(nrow(data)) / 2`

    In the mean time, a numeric value can be passed to `beta` as well to
    explicitly specify the penalty.
