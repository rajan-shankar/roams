# Extract Attributes from a ROAMS SSM List

Retrieves a specified attribute from each model within a
`roams_SSM_list` object. Also works if a single model of class
`roams_SSM` is provided.

## Usage

``` r
get_attribute(model_list, attribute)
```

## Arguments

- model_list:

  An object of class `roams_SSM_list` or a single `roams_SSM` model. May
  optionally include in-sample information added via
  `attach_insample_info`.

- attribute:

  A character string specifying the name of the attribute to extract.
  Must be one of the available scalar attributes (e.g. `BIC`, `lambda`)
  or list/vector attributes (e.g. `filtered_states`, `smoothed_states`).

## Value

If `attribute` is a scalar attribute (e.g., `BIC` or `lambda`), returns
a numeric or character vector containing that attribute across all
models.

If `attribute` is a list-like attribute (e.g., `par` or `gamma`),
returns a list of those values, across all models.

## Details

Available attributes depend on whether in-sample information has been
attached using
[`attach_insample_info()`](https://rajan-shankar.github.io/roams/reference/attach_insample_info.md).
If not, only core model components (e.g. `par`, `gamma`, `y`) and scalar
metrics (e.g. `BIC`, `loglik`) are available.

**Scalar Attributes (always available):**

- `lambda`: Regularization penalty value.

- `prop_outlying`: Proportion of outlying time points identified.

- `BIC`: Bayesian Information Criterion.

- `loglik`: Log-likelihood.

- `RSS`: Residual sum of squares.

- `iterations`: Number of IPOD iterations.

- `value`: Final objective function value.

- `convergence`: Convergence status of optimizer.

- `message`: Optimizer termination message.

**List/Vector Attributes that are always available:**

- `par`, `gamma`, `y`, `counts`

**List/Vector Attributes that are only available if in-sample info is
attached:**

- `smoothed_states`, `filtered_states`, `predicted_states`

- `smoothed_observations`, `filtered_observations`,
  `predicted_observations`

- `smoothed_states_var`, `filtered_states_var`, `predicted_states_var`,
  `predicted_observations_var`

- `mahalanobis_residuals`

Note that these 'in-sample info' attributes are typically only available
if `model_list` is a single `roams_SSM` with in-sample information
already attached using
[`attach_insample_info()`](https://rajan-shankar.github.io/roams/reference/attach_insample_info.md).
If `model_list` is a `roams_SSM_list`, these attributes will not be
available unless all models in the list have in-sample information
attached.

## See also

[`attach_insample_info`](https://rajan-shankar.github.io/roams/reference/attach_insample_info.md),
[`roams_SSM`](https://rajan-shankar.github.io/roams/reference/roams_SSM.md)
