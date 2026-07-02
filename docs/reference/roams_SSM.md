# Robust Outlier-Adjusted Mean-Shift Estimation of State Space Models

Fits a robust state space model to multivariate time series data using
iterative parameter estimation and outlier detection. This procedure is
inspired by the Iterative Procedure for Outlier Detection (IPOD)
algorithm of She and Owen (2011) and is applied over a sequence of
regularization parameters (\\\lambda\\'s), identifying outliers via
Mahalanobis residuals and re-fitting the model iteratively.

## Usage

``` r
roams_SSM(
  y,
  init_par,
  build,
  build_args = list(),
  num_lambdas = 20,
  custom_lambdas = NA,
  cores = 1,
  B = 50,
  lower = NA,
  upper = NA,
  tol = 1e-04,
  lambda_min = 2,
  excessive_outliers_iter_limit = 1,
  control = list(parscale = init_par),
  thresholding = "hard"
)
```

## Arguments

- y:

  A numeric matrix of observations, with each row corresponding to a
  time point.

- init_par:

  A numeric vector of initial parameter values for optimization.

- build:

  A function whose first argument is a parameter vector and that returns
  a `dlm` model (as used in
  [`dlm::dlmMLE()`](https://rdrr.io/pkg/dlm/man/dlmMLE.html)). The
  `specify_SSM` function can be used to create this `build` function.
  Additional named arguments can be declared in `build` and supplied via
  `build_args`.

- build_args:

  An optional named list of additional arguments to forward to `build`
  on every call. For example, `build_args = list(u = covariate_matrix)`
  allows `build` to be written as `function(par, u) { ... }` without
  requiring a closure or factory function. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- num_lambdas:

  Integer. The number of \\\lambda\\ values to evaluate. Ignored if
  `custom_lambdas` is specified. Default is 20.

- custom_lambdas:

  Optional numeric vector. If supplied, these are the exact \\\lambda\\
  values used for model fitting. If not provided or set to `NA`, then
  `num_lambdas` \\\lambda\\'s are automatically chosen.

- cores:

  Integer. Number of CPU cores to use for parallel processing. Default
  is 1 (sequential execution).

- B:

  Integer. Maximum number of iterations per \\\lambda\\. Default is 50.

- lower:

  Optional numeric vector of lower bounds for optimization. If `NA`,
  defaults to `-Inf` for all parameters. Must be of same length as
  `init_par`.

- upper:

  Optional numeric vector of upper bounds for optimization. If `NA`,
  defaults to `Inf` for all parameters. Must be of same length as
  `init_par`.

- tol:

  Tolerance level for checking convergence of the ROAMS procedure.
  Default is `1e-4`.

- lambda_min:

  Minimum \\\lambda\\ value to consider when constructing the sequence
  of \\\lambda\\'s. Ignored if `custom_lambdas` is specified. Default is
  2.

- excessive_outliers_iter_limit:

  Integer. Maximum number of iterations allowed where \\\ge 50\\\\ of
  timepoints are flagged as outliers. This many outliers suggests
  \\\lambda\\ is too low. Allows ROAMS to get through these
  \\\lambda\\'s quicker. Default is 1.

- control:

  A named list of control options to pass to `optim` via
  [`dlm::dlmMLE()`](https://rdrr.io/pkg/dlm/man/dlmMLE.html). Default is
  `list(parscale = init_par)`, which can help the optimizer if
  parameters are on vastly different scales.

- thresholding:

  Character string specifying the outlier thresholding rule. Either
  `"hard"` (default) or `"soft"`. Hard thresholding fully removes
  flagged observations (sets them to missing), corresponding to an L0
  penalty on the outlier shifts \\\gamma\\. Soft thresholding applies a
  continuous LASSO shrinkage to each \\\gamma_t\\, partially adjusting
  observations rather than discarding them entirely, which corresponds
  to an L1 penalty.

## Value

If more than one \\\lambda\\ values are used, returns an object of class
`roams_SSM_list` — a list containing a `roams_SSM` model for each
\\\lambda\\. If only one \\\lambda\\ value is used (i.e.
`custom_lambdas` is manually specified as a single value), returns a
single `roams_SSM` object.

Each `roams_SSM` object includes:

- `lambda` - The \\\lambda\\ value used.

- `prop_outlying` - Proportion of non-missing time points identified as
  outliers.

- `BIC` - Bayesian Information Criterion of the final model.

- `loglik` - Log-likelihood of the fitted model.

- `RSS` - Residual sum of squares.

- `gamma` - Matrix of estimated outlier adjustments.

- `iterations` - Number of iterations performed.

- Optimization output from
  [`dlm::dlmMLE()`](https://rdrr.io/pkg/dlm/man/dlmMLE.html) from the
  final iteration.

- `y` - The original data matrix.

- `build` - The original build function used to specify the model.

## Details

The ROAMS procedure alternates between estimating model parameters via
maximum likelihood and identifying outlying observations based on
Mahalanobis distance of residuals. For each iteration:

1.  A `dlm` model is fit using
    [`dlm::dlmMLE()`](https://rdrr.io/pkg/dlm/man/dlmMLE.html).

2.  Mahalanobis distance of residuals (Mahalanobis residuals) are
    computed.

3.  Outlier shift estimates \\\gamma_t\\ are updated using the selected
    thresholding rule.

4.  The adjusted data `adj_y` is updated and passed to the next
    iteration.

Under **hard thresholding**, observations whose penalised Mahalanobis
score \\d_t^2 + \log\|\mathbf{S}\_{t\|t-1}\|\\ exceeds \\\lambda^2\\ are
fully removed (set to missing), equivalent to setting \\\hat{\gamma}\_t
= r_t\\. Under **soft thresholding**, a group-LASSO shrinkage operator
is applied: \\\hat{\gamma}\_t = \max(1 - \lambda / d_t,\\ 0)\\ r_t\\,
where \\d_t\\ is the Mahalanobis distance of the one-step-ahead
residual. Partial adjustments are retained in the data
(`adj_y = y - gamma`) rather than treating observations as missing.

The algorithm stops when the change in parameters and outlier estimates
is sufficiently small or if too many outliers are detected (more than
50% of complete observations).

## References

She, Y., & Owen, A. B. (2011). Outlier Detection Using Nonconvex
Penalized Regression. *Journal of the American Statistical Association,
106*(494), 626–639. https://doi.org/10.1198/jasa.2011.tm10390

## See also

[`dlmMLE`](https://rdrr.io/pkg/dlm/man/dlmMLE.html),
[`best_BIC_model`](https://rajan-shankar.github.io/roams/reference/best_BIC_model.md),
[`outlier_target_model`](https://rajan-shankar.github.io/roams/reference/outlier_target_model.md),
[`get_attribute`](https://rajan-shankar.github.io/roams/reference/get_attribute.md),
[`autoplot.roams_SSM_list`](https://rajan-shankar.github.io/roams/reference/autoplot.roams_SSM_list.md),
[`attach_insample_info`](https://rajan-shankar.github.io/roams/reference/attach_insample_info.md),
[`oos_filter`](https://rajan-shankar.github.io/roams/reference/oos_filter.md),
[`specify_SSM`](https://rajan-shankar.github.io/roams/reference/specify_SSM.md)
