# Classical State Space Model Fit

Fits a state space model using classical maximum likelihood estimation
with no attempt to detect or account for outliers. This serves as a
baseline model for comparison.

## Usage

``` r
classical_SSM(
  y,
  init_par,
  build,
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

An object of class `classical_SSM` containing the optimization result,
the original data, and the original build function.

## See also

[`dlmMLE`](https://rdrr.io/pkg/dlm/man/dlmMLE.html),
[`roams_SSM`](https://rajan-shankar.github.io/roams/reference/roams_SSM.md),
[`attach_insample_info`](https://rajan-shankar.github.io/roams/reference/attach_insample_info.md),
[`oos_filter`](https://rajan-shankar.github.io/roams/reference/oos_filter.md),
[`specify_SSM`](https://rajan-shankar.github.io/roams/reference/specify_SSM.md)
