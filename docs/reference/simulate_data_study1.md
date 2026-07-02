# Simulate DCRW Data for Study 1: Different Outlier Configurations

Simulates data sets under a first-difference correlated random walk
(DCRW) state-space model for Study 1 in the paper. This study evaluates
performance under different types of outlier configurations. All
arguments have default values matching the simulation setup used in the
paper.

## Usage

``` r
simulate_data_study1(
  sample_sizes = c(100, 200, 500, 1000),
  samples = 100,
  n_oos = 20,
  contamination = 0.1,
  distance = 5,
  sd_cluster = 2,
  mean_cluster = c(20, 20),
  multi_level_distances = c(distance - 2, distance, distance + 2),
  phi_coef = 0.8,
  sigma2_w_lon = 0.1,
  sigma2_w_lat = 0.1,
  sigma2_v_lon = 0.4,
  sigma2_v_lat = 0.4,
  initial_state = c(0, 0, 0, 0),
  seed = NA
)
```

## Arguments

- sample_sizes:

  Vector of sample sizes \\n\\ for each simulated dataset. Default is
  `c(100, 200, 500, 1000)`.

- samples:

  Number of simulated data sets per \\n\\ and per outlier configuration.
  Default is 100.

- n_oos:

  Number of out-of-sample (future) timesteps. Default is 20.

- contamination:

  Proportion of contaminated (outlying) observations. Default is 0.1.

- distance:

  Distance used for fixed-distance outliers. Default is 5.

- sd_cluster:

  Standard deviation of cluster for cluster-based outliers. Default is
  2.

- mean_cluster:

  Mean vector of cluster for cluster-based outliers. Default is
  `c(20, 20)`.

- multi_level_distances:

  Vector of distances for multi-level outliers. Must be of length 3.
  Default is `c(3, 5, 7)`.

- phi_coef:

  Autocorrelation parameter in the DCRW transition matrix. Ranges
  between 0 and 1. Default is 0.8.

- sigma2_w_lon, sigma2_w_lat:

  State noise variances (longitude and latitude). Default is 0.1 each.

- sigma2_v_lon, sigma2_v_lat:

  Observation noise variances (longitude and latitude). Default is 0.4
  each.

- initial_state:

  Initial state vector of length 4. Default is `c(0, 0, 0, 0)`.

- seed:

  Optional random seed for reproducibility. Default is `NA`. Use
  `seed = 1302` to reproduce the same data as in the paper.

## Value

A tibble containing the simulated data sets. Each row corresponds to a
simulated data set and includes fields for sample size, outlier
configuration (setting), outliers, clean data, noisy observations, and
out-of-sample values.

## See also

[`simulate_data_study2`](https://rajan-shankar.github.io/roams/reference/simulate_data_study2.md)

## Examples

``` r
data_study1 = simulate_data_study1(samples = 5, seed = 123)
```
