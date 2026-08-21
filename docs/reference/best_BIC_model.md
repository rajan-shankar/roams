# Select Best Model Based on an Information Criterion

Extracts the best-fitting model from a `roams_SSM_list` object according
to a chosen information criterion, while excluding models with more than
50% outlying observations.

## Usage

``` r
best_BIC_model(model_list, ic = c("BIC", "AIC", "HQIC"), k = NULL)
```

## Arguments

- model_list:

  An object of class `roams_SSM_list`.

- ic:

  Character string specifying the information criterion. One of `"BIC"`
  (default), `"AIC"`, or `"HQIC"` (Hannan-Quinn). Ignored if `k` is
  supplied.

- k:

  Optional numeric penalty multiplier per model parameter, as in
  [`step`](https://rdrr.io/r/stats/step.html). Overrides `ic` when
  non-`NULL`. Common values: `2` (AIC), `log(n)` (BIC),
  `2 * log(log(n))` (HQIC).

## Value

A single `roams_SSM` object corresponding to the model with the smallest
value of the chosen information criterion among those with fewer than
50% outlying observations.

## Details

The criterion is computed as \\-2\ell + k \cdot p\\, where \\\ell\\ is
the log-likelihood, \\p\\ is the number of non-zero components, and
\\k\\ is the penalty multiplier. The `ic` argument provides named
shortcuts; `k` allows a fully custom penalty and overrides `ic` when
non-`NULL`.

## Examples

``` r
if (FALSE) { # \dontrun{
# BIC (default)
best_BIC_model(models)

# Named shortcuts
best_BIC_model(models, ic = "AIC")
best_BIC_model(models, ic = "HQIC")

# Custom penalty via k (overrides ic)
best_BIC_model(models, k = 2)
} # }
```
