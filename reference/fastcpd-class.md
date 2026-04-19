# An S4 class to store the output created with [`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md)

This S4 class stores the output from
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) and
[fastcpd.family](https://fastcpd.xingchi.li/reference/fastcpd_family.md).
A fastcpd object consist of several slots including the call to
[`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md), the data
used, the family of the model, the change points, the cost values, the
residuals, the estimated parameters and a boolean indicating whether the
model was fitted with only change points or with change points and
parameters, which you can select using `@`.

## Slots

- `call`:

  The call of the function.

- `data`:

  The data passed to the function.

- `order`:

  The order of the time series model.

- `family`:

  The family of the model.

- `cp_set`:

  The set of change points.

- `cost_values`:

  The cost function values for each segment.

- `residuals`:

  The residuals of the model with change points. Used only for built-in
  families.

- `thetas`:

  The estimated parameters for each segment. Used only for built-in
  families.

- `cp_only`:

  A boolean indicating whether
  [`fastcpd()`](https://fastcpd.xingchi.li/reference/fastcpd.md) was run
  to return only the change points or the change points with the
  estimated parameters and cost values for each segment.
