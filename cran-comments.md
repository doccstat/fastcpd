## Submission notes

The DOI in the CITATION is for a JSS publication that will be registered
after publication on CRAN.

## Updates since last CRAN release (0.16.2)

### fastcpd 1.0.0

*   Add `exponential` family (`fastcpd_exponential` / `fastcpd.exponential`)
    for piecewise-constant rate change detection in exponentially distributed
    data.
*   Add support for compiled C++ cost functions via `Rcpp::XPtr` in
    `family = "custom"`, eliminating R-call overhead in the hot PELT/SeGD
    loop. All three parameters (`cost`, `cost_gradient`, `cost_hessian`) accept
    tagged external pointers alongside R closures; arity is conveyed via
    `attr(xptr, "fastcpd_cost_arity")`.
*   Remove bundled third-party source files (`ref_fastglm_*`, `ref_tseries*`);
    GLM and GARCH fitting now uses internal implementations (`FitGlm`,
    `FitGarch`) with no external source copies.
*   Move `fastcpd_impl.h` from `inst/include/` to `src/`; remove
    `inst/include/` directory entirely.

### fastcpd 0.99.9

*   Add JSS publication reference: Li and Zhang (2026)
    <doi:10.18637/jss.v116.i06>.
*   Update CITATION to include the Journal of Statistical Software article.

### fastcpd 0.20.0

*   Use templates generation.

### fastcpd 0.19.0

*   Relicense from GPL-3.0 to Apache 2.0.

### fastcpd 0.18.0

*   Release wheel for Python on Linux and macOS.

### fastcpd 0.17.0

*   Introduce Python version.

### fastcpd 0.16.3

*   Move `fastcpd` to the top in the README.
