# Compute Out-of-Sample Inference for Fitted State Space Model

Applies the fitted model parameters to a user-supplied out-of-sample
dataset to compute predicted and filtered states and observations.
Robust and classical inference procedures are supported depending on the
class of the input model.

## Usage

``` r
oos_filter(
  y_oos,
  model,
  build,
  build_args = list(),
  outlier_locs = rep(0, nrow(y_oos)),
  threshold = sqrt(qchisq(0.99, ncol(y_oos))),
  multiplier = 2
)
```

## Arguments

- y_oos:

  A numeric matrix containing out-of-sample observations. Each row
  corresponds to a time point.

- model:

  A fitted model object of class `roams_SSM`, `classical_SSM`,
  `oracle_SSM`, `huber_robust_SSM`, or `trimmed_robust_SSM`.

- build:

  A function whose first argument is a parameter vector and that returns
  a `dlm` model object. The `specify_SSM` function can be used to create
  this `build` function.

- build_args:

  An optional named list of additional arguments to forward to `build`.
  Should match what was supplied to the fitting function (e.g.
  `roams_SSM`) unless a different covariate matrix is needed for the
  out-of-sample period. Default is
  [`list()`](https://rdrr.io/r/base/list.html).

- outlier_locs:

  A logical or binary vector of the same length as `nrow(y)`, indicating
  time points to be treated as missing (i.e., time points that are known
  to be outliers). Used only with `oracle_SSM` models.

- threshold:

  Mahalanobis distance threshold for detecting out-of-sample outliers in
  `roams_SSM` models. Set to `Inf` to recover the usual Kalman filter.
  Default is `sqrt(qchisq(0.99, ncol(y)))`.

- multiplier:

  Multiplier for how quickly the filter grows its filtered state
  variance (uncertainty) after detecting an outlier in `roams_SSM`
  models. It is the tuneable parameter \\b\\ of the \`fast-updating
  threshold' filter. Only works if `threshold` is not `Inf`. Default is
  `2`.

## Value

A named list containing out-of-sample inference results:

- `filtered_states`:

  Filtered state estimates using out-of-sample data.

- `predicted_states`:

  One-step-ahead state predictions.

- `filtered_observations`:

  Expected observations given past out-of-sample data.

- `predicted_observations`:

  One-step-ahead forecasts of observations.

- `filtered_states_var`:

  List of filtered state variance matrices.

- `predicted_states_var`:

  List of one-step-ahead state prediction variances.

- `predicted_observations_var`:

  List of one-step-ahead observation forecast variances.

- `mahalanobis_residuals`:

  Vector of Mahalanobis distances of residuals from predicted
  observations.

- `outliers_flagged`:

  Vector of 1's and 0's indicating whether timepoints are flagged as
  outlying or not based on the `threshold` supplied (only available if
  `model` is of class `roams_SSM`).

## Details

The function reuses the model's fitted parameters to generate inference
on new data `y_oos`. Robust variants use appropriate robust filters,
while the classical and oracle models use standard Kalman filtering. For
`oracle_SSM` models, observations flagged in `outlier_locs` are treated
as missing during filtering.

## See also

[`attach_insample_info`](https://rajan-shankar.github.io/roams/reference/attach_insample_info.md),
[`specify_SSM`](https://rajan-shankar.github.io/roams/reference/specify_SSM.md)
