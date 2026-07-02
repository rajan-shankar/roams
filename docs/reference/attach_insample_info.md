# Attach In-Sample Information to Fitted State Space Model

Attaches detailed in-sample information—such as predicted, filtered, and
smoothed states and observations—to a model object fitted using any of
the package’s supported SSM estimation methods. These quantities are not
stored by default in model objects due to their potentially large memory
footprint.

## Usage

``` r
attach_insample_info(model)
```

## Arguments

- model:

  A fitted model object of class `roams_SSM`, `classical_SSM`,
  `oracle_SSM`, `huber_robust_SSM`, or `trimmed_robust_SSM`.

## Value

A modified version of the input model object, with an additional class
`insample_info`, and the following in-sample elements appended:

- `filtered_states`:

  Filtered state estimates using data up to each time point.

- `predicted_states`:

  One-step-ahead state predictions.

- `filtered_observations`:

  Expected observations given data up to each time point.

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

For models of class `roams_SSM`, `classical_SSM`, or `oracle_SSM`, the
following additional elements are also attached:

- `smoothed_states`:

  Posterior means of hidden states using all data.

- `smoothed_observations`:

  Posterior mean of the observed series based on smoothed states.

- `smoothed_states_var`:

  List of smoothed state variance matrices.

These smoothed attributes are obtained using the RTS smoothing algorithm
(Rauch et al. 1965).

## Details

The attached outputs enable richer diagnostics, outlier inspection, and
plotting. For `huber_robust_SSM` and `trimmed_robust_SSM` models,
in-sample information is computed using a custom robust filtering
function, and smoothed quantities (`smoothed_states`,
`smoothed_observations`, and `smoothed_states_var`) are **not
available**. This function should only be applied once to a model
object.

## References

Rauch, H.E., Tung, F., Striebel, C.T. (1965). Maximum likelihood
estimates of linear dynamic systems. *AIAA Journal 3*(8), 1445–1450.
https://doi.org/10.2514/3.3166

## See also

[`oos_filter`](https://rajan-shankar.github.io/roams/reference/oos_filter.md)
