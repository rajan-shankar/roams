# Autoplot for a ROAMS SSM List

Generates diagnostic plots for a list of robust state space models fit
across a sequence of \\\lambda\\ values. Two model attributes are
plotted against \\\lambda\\, each in its own panel, with a vertical
dashed line indicating the model with the lowest BIC among those with
fewer than 50% outliers.

## Usage

``` r
# S3 method for class 'roams_SSM_list'
autoplot(object, attribute1 = "BIC", attribute2 = "prop_outlying", ...)
```

## Arguments

- object:

  An object of class `roams_SSM_list` as returned by
  [`roams_SSM`](https://rajan-shankar.github.io/roams/reference/roams_SSM.md)
  when multiple \\\lambda\\ values are used.

- attribute1:

  A character string indicating the first model attribute to plot on the
  top panel. Options include `"lambda"`, `"prop_outlying"`, `"BIC"`,
  `"AIC"`, `"HQIC"`, `"loglik"`, `"RSS"`, `"iterations"`, and `"value"`.
  Defaults to `"BIC"`.

- attribute2:

  A character string indicating the second model attribute to plot on
  the bottom panel. Uses the same options as `attribute1`. Defaults to
  `"prop_outlying"`.

- ...:

  Other arguments passed to specific methods. Not used in this method.

## Value

A `ggplot` object arranged using the patchwork package, showing the
specified attributes plotted against \\\lambda\\ in two vertically
stacked panels.

## Details

In each panel, a red dashed vertical line indicates the model with the
lowest BIC among those with fewer than 50% outlying time points, serving
as a heuristic for robust model selection.

## See also

[`roams_SSM`](https://rajan-shankar.github.io/roams/reference/roams_SSM.md),
[`get_attribute`](https://rajan-shankar.github.io/roams/reference/get_attribute.md)
