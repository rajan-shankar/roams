# Select Best Model Based on BIC

Extracts the best-fitting model from a `roams_SSM_list` object according
to the Bayesian Information Criterion (BIC), while excluding models with
more than 50% outlying observations.

## Usage

``` r
best_BIC_model(model_list)
```

## Arguments

- model_list:

  An object of class `roams_SSM_list`

## Value

A single `roams_SSM` object corresponding to the model with the smallest
BIC among those with fewer than 50% outlying observations.

## Examples

``` r
if (FALSE) { # \dontrun{
# Assuming `models` is a roams_SSM_list:
best_model <- best_BIC_model(models)
} # }

```
