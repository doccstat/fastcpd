## Updates since last CRAN release (0.10.3)

### fastcpd 0.11.1

*   Add package comparison with `gfpop`, `jointseg` and `Rbeast`.

### fastcpd 0.11.0

*   **Note**: From now on, MBIC is used as the default penalty selection for
    `beta` parameter.
*   Add penalty selection criteria using

    1. BIC: `(p + 1) * log(nrow(data)) / 2`
    1. Modified BIC: `(p + 2) * log(nrow(data)) / 2` with adjusted cost
       function.
    1. MDL: `(p + 2) * log(nrow(data)) / 2` with adjusted cost function.

    In the mean time, a numeric value can be passed to `beta` as well to
    explicitly specify the penalty for BIC.
