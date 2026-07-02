# roams

**roams** provides functions for estimating parameters and detecting
outliers in state-space models using the ROAMS (Robust Outlier-Adjusted
Mean-Shift) framework. It includes functionality for fitting benchmark
models, visualizing model selection criteria, and evaluating in-sample
and out-of-sample performance. Simulation tools for generating synthetic
data under a first-differenced correlated random walk (DCRW) model are
also included. Designed with flexibility for user-defined model
structures via the
[`specify_SSM()`](https://rajan-shankar.github.io/roams/reference/specify_SSM.md)
interface.

See the corresponding ROAMS paper (pre-print) here:
[arxiv.org/abs/2511.15155](https://arxiv.org/abs/2511.15155)

## Installation

You can install the development version of **roams** from GitHub with:

``` r

# install.packages("pak")
pak::pak("rajan-shankar/roams")
```

OR if using Windows:

``` r

# install.packages("devtools")
devtools::install_github("rajan-shankar/roams")
```
