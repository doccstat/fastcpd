## Updates since last CRAN release (0.10.3)

### fastcpd 0.12.1

*   Remove useless C++ codes.
*   Add more debug points in C++.
*   Add more examples for data `well_log`.
*   Add detection comparison for `well_log` data.
*   Add a variance estimator for median change.
*   Deprecate `winsorize_minval` and `winsorize_maxval`.

### fastcpd 0.12.0

*   Add Rice estimation for ARMA model variance estimation.
*   Add time comparison using Well-log data in vignettes.

### fastcpd 0.11.3

*   Add Rice estimator for mean change variance estimation.

### fastcpd 0.11.2

*   Export variance estimator function for linear models.

### fastcpd 0.11.1

*   Add package comparison with `CptNonPar`, `gfpop`, `InspectChangepoint`,
    `jointseg`, `Rbeast` and `VARDetect`.

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
