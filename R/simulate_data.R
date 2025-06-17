#' Simulate DCRW Data for Study 1: Different Outlier Configurations
#'
#' Simulates datasets under a first-difference correlated random walk (DCRW) state-space model
#' for Study 1 in the paper. This study evaluates performance under different types of outlier contamination.
#' All arguments have default values matching the simulation setup used in the paper.
#'
#' @param sample_sizes Vector of sample sizes for each simulated dataset. Default is \code{c(100, 200, 500)}.
#' @param samples Number of Monte Carlo replicates per scenario. Default is 100.
#' @param n_oos Number of out-of-sample (future) timesteps. Default is 100.
#' @param contamination Proportion of contaminated (outlying) observations. Default is 0.1.
#' @param distance Distance used for fixed-distance outliers. Default is 5.
#' @param sd_cluster Standard deviation of cluster-based outliers. Default is 2.
#' @param mean_cluster Mean vector for cluster-based outliers. Default is \code{c(20, 20)}.
#' @param multi_level_distances Vector of distances for multi-level outliers. Must be of length 3. Default is \code{c(3, 5, 7)}.
#' @param phi_coef AR(1) parameter in the DCRW transition matrix. Default is 0.8.
#' @param sigma2_w_lon, sigma2_w_lat Process noise variances (longitude and latitude). Default is 0.1 each.
#' @param sigma2_v_lon, sigma2_v_lat Observation noise variances (longitude and latitude). Default is 0.4 each.
#' @param initial_state Initial state vector of length 4. Default is \code{c(0, 0, 0, 0)}.
#' @param seed Optional random seed for reproducibility. Default is \code{NA}.
#'
#' @return A tibble containing the simulated datasets. Each row corresponds to a simulation replicate and includes
#' fields for sample size, contamination type, outliers, clean data, noisy observations, and out-of-sample values.
#'
#' @examples
#' data_study1 <- simulate_data_study1(samples = 10, seed = 123)
#'
#' @export
simulate_data_study1 = function(
  sample_sizes = c(100, 200, 500),
  samples = 100,
  n_oos = 100,
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
  seed = NA  # set to 1302 for same data as paper
  ) {
  # This function simulates data for the first simulation study in the paper "Study 1: Different outlier configurations".

  # If incorrect argument lengths, then stop the function with an error:
  if (length(n_oos) > 1 || length(contamination) > 1 ||
      length(distance) > 1 || length(sd_cluster) > 1 || length(mean_cluster) != 2 ||
      length(multi_level_distances) != 3 || length(phi_coef) > 1 || length(sigma2_w_lon) > 1 ||
      length(sigma2_w_lat) > 1 || length(sigma2_v_lon) > 1 || length(sigma2_v_lat) > 1 ||
      length(initial_state) != 4) {
    stop("Check the lengths of your arguments. Argument 'mean_cluster' must be of length 2, argument 'multi_level_distances' must be of length 3, argument 'initial_state' must be of length 4, and argument 'sample_sizes' must be of length >= 1. All other arguments must be of length 1.")
  }

  samples = 1:samples
  ns = sample_sizes
  multi_level_dists = multi_level_distances
  multi_level_dists = multi_level_dists[order(multi_level_dists)]
  dist = distance
  cont = contamination

  if (!is.na(seed)) {
    set.seed(seed)
  }

  data_sets = tibble(
    sample = numeric(),
    n = numeric(),
    setting = character(),
    y_oos = list(),
    x_oos = list(),
    y = list(),
    y_clean = list(),
    x = list(),
    outlier_locs = list(),
    outlier_levels = list()
  )

  settings = c("clean data", "fixed distance", "multi-level", "cluster")

  # DCRW data-generating process
  Phi = diag(c(1+phi_coef, 1+phi_coef, 0, 0))
  Phi[1,3] = -phi_coef
  Phi[2,4] = -phi_coef
  Phi[3,1] = 1
  Phi[4,2] = 1

  A = diag(4)[1:2,]
  Q = diag(c(sigma2_w_lon, sigma2_w_lat, 0, 0))
  R = diag(c(sigma2_v_lon, sigma2_v_lat))

  x0 = initial_state

  # Simulate data loop
  for (setting in settings) {
    for (n in ns) {
      for (sample_number in samples) {
        # State process
        x = matrix(nrow = 4, ncol = n + n_oos)
        x[,1] = Phi %*% x0 + MASS::mvrnorm(mu = rep(0,4), Sigma = Q)
        for (t in 2:(n + n_oos)) {
          x[,t] = Phi %*% x[,t-1] + MASS::mvrnorm(mu = rep(0,4), Sigma = Q)
        }

        # Observation process
        y_clean = matrix(nrow = 2, ncol = n + n_oos)
        for (t in 1:(n + n_oos)) {
          y_clean[,t] = A %*% x[,t] + MASS::mvrnorm(mu = rep(0,2), Sigma = R)
        }

        # Generate outliers
        generate_outliers = function(setting, n, cont, dist, multi_level_dists, sd_cluster) {
          outlier_levels = NA

          if (setting == "clean data") {
            outlier_locs = rep(0, n)
            outlier_vecs = matrix(0, nrow = 2, ncol = n)

          } else {
            outlier_locs = c(
              0,  # first point cannot be an outlier
              sample(c(
                rep(1, round(n*cont)),
                rep(0, n - round(n*cont) - 1)
              ))
            )

            angles = runif(n, max = 2*pi)
            outlier_unit_vecs = do.call(
              cbind,
              map(angles, ~ c(cos(.), sin(.)))
            )

            if (setting == "fixed distance") {
              outlier_vecs = outlier_unit_vecs*dist

            } else if (setting == "multi-level") {
              levels = sample(c(
                rep(multi_level_dists[1], floor(n*cont / 3)),
                rep(multi_level_dists[2], floor(n*cont / 3)),
                rep(multi_level_dists[3], n*cont - 2*floor(n*cont / 3))
              ))

              counter = 1
              for (t in 1:n) {
                if (outlier_locs[t] == 1) {
                  outlier_levels[t] = levels[counter]
                  counter = counter + 1
                } else {
                  outlier_levels[t] = 0
                }
              }

              outlier_vecs = outlier_unit_vecs %*% diag(outlier_levels)

            } else if (setting == "cluster") {
              outlier_vecs = matrix(0, nrow = 2, ncol = n)
            }
          }

          return(list(
            "locs" = outlier_locs,
            "levels" = outlier_levels,
            "vecs" = outlier_vecs
          ))
        }
        outlier_info = generate_outliers(setting, n, cont, dist, multi_level_dists, sd_cluster)

        # Observational outlier process
        y = y_clean
        if (setting != "cluster") {
          for (t in 1:n) {
            y[,t] = y_clean[,t] + outlier_info$locs[t]*outlier_info$vecs[,t]
          }
        } else {
          for (t in 1:n) {
            y[,t] = (1 - outlier_info$locs[t])*y_clean[,t] + outlier_info$locs[t]*(MASS::mvrnorm(mu = mean_cluster, Sigma = sd_cluster^2*diag(2)))
          }
        }

        data_sets = data_sets %>%
          add_row(
            sample = sample_number,
            n = n,
            setting = setting,
            y_oos = list(y[,(n + 1):(n + n_oos)]),
            x_oos = list(x[,(n + 1):(n + n_oos)]),
            y = list(y[,1:n]),
            y_clean = list(y_clean[,1:n]),
            x = list(x[,1:n]),
            outlier_locs = list(outlier_info$locs),
            outlier_levels = list(outlier_info$levels)
          )

      }
    }
  }

  data_sets = data_sets %>%
    select(sample, n, setting, y, x, y_oos, x_oos, y_clean, outlier_locs, outlier_levels)

  return(data_sets)
}

#' Simulate DCRW Data for Study 2: Increasing Contamination and Varying Outlier Distance
#'
#' Simulates datasets under a first-difference correlated random walk (DCRW) state-space model
#' for Study 2 in the paper. This study examines the impact of increasing contamination levels
#' and varying outlier distances on model performance. All arguments have default values matching
#' the simulation setup used in the paper.
#'
#' @param samples Number of Monte Carlo replicates per scenario. Default is 100.
#' @param n Number of in-sample timesteps. Default is 200.
#' @param max_contamination Maximum proportion of contaminated (outlying) observations. Default is 0.2.
#' @param distances Vector of five distances for additive outliers. Must be of length 5. Default is \code{c(1, 3, 5, 7, 9)}.
#' @param n_oos Number of out-of-sample (future) timesteps. Default is 100.
#' @param phi_coef AR(1) parameter in the DCRW transition matrix. Default is 0.8.
#' @param sigma2_w_lon, sigma2_w_lat Process noise variances (longitude and latitude). Default is 0.1 each.
#' @param sigma2_v_lon, sigma2_v_lat Observation noise variances (longitude and latitude). Default is 0.4 each.
#' @param initial_state Initial state vector of length 4. Default is \code{c(0, 0, 0, 0)}.
#' @param seed Optional random seed for reproducibility. Default is \code{NA}.
#'
#' @return A tibble containing the simulated datasets. Each row corresponds to a simulation replicate and includes
#' fields for contamination level, distance, outliers, clean data, noisy observations, and out-of-sample values.
#'
#' @examples
#' data_study2 <- simulate_data_study2(samples = 10, seed = 456)
#'
#' @export
simulate_data_study2 = function(
  samples = 100,
  n = 200,
  max_contamination = 0.2,
  distances = c(1, 3, 5, 7, 9),
  n_oos = 100,
  phi_coef = 0.8,
  sigma2_w_lon = 0.1,
  sigma2_w_lat = 0.1,
  sigma2_v_lon = 0.4,
  sigma2_v_lat = 0.4,
  initial_state = c(0, 0, 0, 0),
  seed = NA  # set to 205 for same data as paper
  ) {
  # This function simulates data for the second simulation study in the paper "Study 2: Increasing contamination rate and decreasing distance".

  # If incorrect argument lengths, then stop the function with an error:
  if (length(n) > 1 || length(max_contamination) > 1 ||
      length(n_oos) > 1 || length(phi_coef) > 1 || length(sigma2_w_lon) > 1 ||
      length(sigma2_w_lat) > 1 || length(sigma2_v_lon) > 1 || length(sigma2_v_lat) > 1 ||
      length(initial_state) != 4 || length(distances) != 5) {
    stop("Check the lengths of your arguments. Argument 'distances' must be of length 5 and argument 'initial_state' must be of length 4. All other arguments must be of length 1.")
  }

  if (!is.na(seed)) {
    set.seed(seed)
  }

  data_sets = tibble(
    sample = numeric(),
    contamination = numeric(),
    distance = numeric(),
    y_oos = list(),
    x_oos = list(),
    y = list(),
    x = list(),
    y_clean = list(),
    outlier_locs = list())

  contaminations = seq(0, max_contamination, length.out = 5)
  distances = distances[order(distances)]

  # DCRW data-generating process
  Phi = diag(c(1+phi_coef, 1+phi_coef, 0, 0))
  Phi[1,3] = -phi_coef
  Phi[2,4] = -phi_coef
  Phi[3,1] = 1
  Phi[4,2] = 1

  A = diag(4)[1:2,]
  Q = diag(c(sigma2_w_lon, sigma2_w_lat, 0, 0))
  R = diag(c(sigma2_v_lon, sigma2_v_lat))

  x0 = initial_state

  for (i in 1:samples) {

    # State process
    x = matrix(0, nrow = 4, ncol = n+n_oos)
    x[,1] = x0
    for (t in 2:(n+n_oos)) {
      x[,t] = Phi %*% x[,t-1] + MASS::mvrnorm(mu = c(0,0,0,0), Sigma = Q)
    }

    y_noises = t(MASS::mvrnorm(n+n_oos, mu = c(0,0), Sigma = R))
    outlier_locs_all = sample(c(
      rep(2, ceiling(max_contamination * n / 4)),
      rep(3, ceiling(max_contamination * n / 4)),
      rep(4, floor(max_contamination * n / 4)),
      rep(5, floor(max_contamination * n / 4)),
      rep(0, n - ceiling(max_contamination * n / 4)*2 - floor(max_contamination * n / 4)*2)
    ), size = n)

    angles = runif(n, max = 2*pi)

    # Observation process
    y = matrix(0, nrow = 2, ncol = n+n_oos)
    for (t in 1:(n+n_oos)) {
      y[,t] = A %*% x[,t] + y_noises[,t]
    }

    y_clean = y

    for (counter in 2:6) {
      outlier_locs = outlier_locs_all
      outlier_locs[which(!(outlier_locs %in% (counter:6)))] = 0
      outlier_locs[which(outlier_locs %in% (counter:6))] = 1

      if (counter == 4) {
        for (distance in distances) {
          # Additive observational outlier process
          z = matrix(0, nrow = 2, ncol = n)
          for (t in 1:n) {
            z[,t] = y[,t] + outlier_locs[t]*c(cos(angles[t]),
                                              sin(angles[t]))*distance
          }
          data_sets = data_sets %>%
            add_row(
              sample = i,
              contamination = contaminations[6+1-counter],
              distance = distance,
              y_oos = list(y[, (n + 1):(n + n_oos)]),
              x_oos = list(x[, (n + 1):(n + n_oos)]),
              y = list(z[, 1:n]),
              x = list(x[, 1:n]),
              y_clean = list(y_clean[, 1:n]),
              outlier_locs = list(outlier_locs)
            )
        }

      } else {
        distance = distances[median(1:length(distances))]
        # Additive observational outlier process
        z = matrix(0, nrow = 2, ncol = n)
        for (t in 1:n) {
          z[,t] = y[,t] + outlier_locs[t]*c(cos(angles[t]),
                                            sin(angles[t]))*distance
        }
        data_sets = data_sets %>%
          add_row(
            sample = i,
            contamination = contaminations[6+1-counter],
            distance = distance,
            y_oos = list(y[, (n + 1):(n + n_oos)]),
            x_oos = list(x[, (n + 1):(n + n_oos)]),
            y = list(z[, 1:n]),
            x = list(x[, 1:n]),
            y_clean = list(y_clean[, 1:n]),
            outlier_locs = list(outlier_locs)
          )
      }

    }

  }

  data_sets = data_sets %>%
    select(sample, contamination, distance, y, x, y_oos, x_oos, y_clean, outlier_locs)

  return(data_sets)

}

