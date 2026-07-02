# Select Model Based on Target Outlier Proportion

Extracts the model from a `roams_SSM_list` object whose estimated
outlier proportion is closest to a user-specified target.

## Usage

``` r
outlier_target_model(model_list, target)
```

## Arguments

- model_list:

  An object of class `roams_SSM_list`.

- target:

  A numeric value between 0 and 1 indicating the desired proportion of
  outlying observations.

## Value

A single `roams_SSM` object whose estimated outlier proportion is
closest to `target`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Select the model with an outlier proportion closest to 10%
target_model <- outlier_target_model(models, target = 0.1)
} # }
```
