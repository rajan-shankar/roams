# Trimmed-Robust State Space Model Fit

Fits a robust state space model by minimizing a trimmed loss function as
per Crevits and Croux (2018). A fixed proportion of the largest
residuals are excluded from the objective, providing robustness against
extreme outliers. The predicted observations used in the trimmed loss
are computed using the Huber robust filter from Cipra and Romera (1997).

## Usage

``` r
trimmed_robust_SSM(
  y,
  init_par,
  build,
  alpha,
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

- alpha:

  Numeric value in the interval \[0, 1) indicating the trimming
  proportion (i.e., the proportion of data to exclude as outliers).

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

An object of class `trimmed_robust_SSM` containing the optimization
result, trimming level \\\alpha\\, the original data, and the original
build function.

## References

Crevits R. and Croux C. (2018). Robust Estimation of Linear State Space
Models. \*Communications in Statistics: Simulation and Computation\*

Cipra, T., Romera, R. (1997). Kalman filter with outliers and missing
observations. \*Test\* 6, 379–395. https://doi.org/10.1007/BF02564705

## See also

[`huber_robust_SSM`](https://rajan-shankar.github.io/roams/reference/huber_robust_SSM.md),
[`roams_SSM`](https://rajan-shankar.github.io/roams/reference/roams_SSM.md),
[`optim`](https://rdrr.io/r/stats/optim.html),
[`attach_insample_info`](https://rajan-shankar.github.io/roams/reference/attach_insample_info.md),
[`oos_filter`](https://rajan-shankar.github.io/roams/reference/oos_filter.md),
[`specify_SSM`](https://rajan-shankar.github.io/roams/reference/specify_SSM.md)
