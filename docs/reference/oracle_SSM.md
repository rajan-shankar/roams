# Oracle State Space Model Fit

Fits a state space model by treating a known set of outliers as missing
data. This benchmark model assumes prior knowledge of outlier locations
and is intended for comparison with automatic outlier detection
procedures.

## Usage

``` r
oracle_SSM(
  y,
  init_par,
  build,
  outlier_locs,
  build_args = list(),
  lower = NA,
  upper = NA,
  control = list(parscale = init_par)
)
```

## Arguments

- y:

  A numeric matrix of observations (time points in rows).

- init_par:

  A numeric vector of initial parameter values.

- build:

  A function that returns a `dlm` model given a parameter vector. The
  [`specify_SSM()`](https://rajan-shankar.github.io/roams/reference/specify_SSM.md)
  function can be used to create this `build` function.

- outlier_locs:

  An integer or logical vector of length equal to the number of time
  points, indicating locations of known outliers.

- lower:

  Optional numeric vector of lower bounds for parameter estimation.
  Defaults to `-Inf`. Must be of same length as `init_par`.

- upper:

  Optional numeric vector of upper bounds for parameter estimation.
  Defaults to `Inf`. Must be of same length as `init_par`.

- control:

  Optional list of control parameters passed to `optim` via
  [`dlm::dlmMLE()`](https://rdrr.io/pkg/dlm/man/dlmMLE.html). Default is
  `list(parscale = init_par)`, which can help the optimizer if
  parameters are on vastly different scales.

## Value

An object of class `oracle_SSM` containing the optimization result, the
provided outlier locations, the original data, and the original build
function.

## See also

[`dlmMLE`](https://rdrr.io/pkg/dlm/man/dlmMLE.html),
[`roams_SSM`](https://rajan-shankar.github.io/roams/reference/roams_SSM.md),
[`attach_insample_info`](https://rajan-shankar.github.io/roams/reference/attach_insample_info.md),
[`oos_filter`](https://rajan-shankar.github.io/roams/reference/oos_filter.md),
[`specify_SSM`](https://rajan-shankar.github.io/roams/reference/specify_SSM.md)
