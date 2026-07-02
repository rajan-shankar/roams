# Specify a State-Space Model in DLM Format

A helper function for users to construct a state-space model in the
format expected by the `dlm` package, which is the package used
internally for some model fitting functions. This function returns a
named list of model components which can be directly used in a
user-defined `build` function passed to any modeling function in this
package.

## Usage

``` r
specify_SSM(
  state_transition_matrix,
  state_noise_var,
  observation_matrix,
  observation_noise_var,
  init_state_mean,
  init_state_var,
  covariates = NULL,
  obs_covariate_index = NULL,
  state_covariate_index = NULL
)
```

## Arguments

- state_transition_matrix:

  A square matrix specifying the state transition dynamics (`GG`).

- state_noise_var:

  A square matrix specifying the variance of the state noise (`W`).

- observation_matrix:

  A matrix mapping the state to the observations (`FF`).

- observation_noise_var:

  A square matrix specifying the variance of the observation noise
  (`V`).

- init_state_mean:

  A vector specifying the initial mean of the state (`m0`).

- init_state_var:

  A square matrix specifying the initial state covariance (`C0`).

- covariates:

  Optional numeric matrix of exogenous inputs (`X`), with \\n\\ rows
  (one per time point) and \\k\\ columns (one per covariate). Default is
  `NULL`.

- obs_covariate_index:

  Optional integer matrix of the same dimensions as `observation_matrix`
  (`JFF`). A nonzero entry `[i, j] = k` means the \\(i,j)\\ element of
  the observation matrix at time \\t\\ is replaced by
  `covariates[t, k]`. Default is `NULL`.

- state_covariate_index:

  Optional integer matrix of the same dimensions as
  `state_transition_matrix` (`JGG`). A nonzero entry `[i, j] = k` means
  the \\(i,j)\\ element of the state transition matrix at time \\t\\ is
  replaced by `covariates[t, k]`. Default is `NULL`.

## Value

A named list with elements `GG`, `W`, `FF`, `V`, `m0`, and `C0`, and
optionally `X`, `JFF`, and/or `JGG` when covariate arguments are
supplied. Suitable for use in a custom `build` function for modeling or
for online filtering (e.g., using
[`oos_filter`](https://rajan-shankar.github.io/roams/reference/oos_filter.md)).

## Details

The letters in the parentheses in the Arguments section correspond to
the naming convention used in the `dlm` package.

Exogenous inputs (covariates) enter the model through the `X`, `JFF`,
and/or `JGG` components. The matrix `X` (supplied via `covariates`) must
be closed over inside the `build` function so that it is available when
the build function is called during fitting.

## See also

[`roams_SSM`](https://rajan-shankar.github.io/roams/reference/roams_SSM.md),
[`oos_filter`](https://rajan-shankar.github.io/roams/reference/oos_filter.md)

## Examples

``` r
build_function = function(par) {
  phi_coef = par[1]
  Phi = diag(c(1+phi_coef, 1+phi_coef, 0, 0))
  Phi[1,3] = -phi_coef
  Phi[2,4] = -phi_coef
  Phi[3,1] = 1
  Phi[4,2] = 1

  A = diag(4)[1:2,]
  Sigma_W = diag(c(par[2], par[3], 0, 0))
  Sigma_V = diag(c(par[4], par[5]))

  mu0 = rep(0, 4)
  P0 = diag(rep(0, 4))

  specify_SSM(
    state_transition_matrix = Phi,
    state_noise_var = Sigma_W,
    observation_matrix = A,
    observation_noise_var = Sigma_V,
    init_state_mean = mu0,
    init_state_var = P0
  )
}
```
