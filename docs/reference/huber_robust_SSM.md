# Huber-Robust State Space Model Fit

Fits a robust state space model by minimizing a Huber loss objective as
per Crevits and Croux (2018), providing protection against moderate
outliers. The predicted observations used in the Huber loss are computed
using the Huber robust filter from Cipra and Romera (1997).

## Usage

``` r
huber_robust_SSM(
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

  Optional list of control parameters passed to `optim`. Default is
  `list(parscale = init_par)`, which can help the optimizer if
  parameters are on vastly different scales.

## Value

An object of class `huber_robust_SSM` containing the optimization
result, the original data, and the original build function.

## References

Crevits R. and Croux C. (2018). Robust Estimation of Linear State Space
Models. \*Communications in Statistics: Simulation and Computation\*

Cipra, T., Romera, R. (1997). Kalman filter with outliers and missing
observations. \*Test\* 6, 379–395. https://doi.org/10.1007/BF02564705

## See also

[`trimmed_robust_SSM`](https://rajan-shankar.github.io/roams/reference/trimmed_robust_SSM.md),
[`roams_SSM`](https://rajan-shankar.github.io/roams/reference/roams_SSM.md),
[`optim`](https://rdrr.io/r/stats/optim.html),
[`attach_insample_info`](https://rajan-shankar.github.io/roams/reference/attach_insample_info.md),
[`oos_filter`](https://rajan-shankar.github.io/roams/reference/oos_filter.md),
[`specify_SSM`](https://rajan-shankar.github.io/roams/reference/specify_SSM.md)
