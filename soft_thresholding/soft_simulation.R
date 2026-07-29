num_cores = 25
num_samples = 100
data_path = "/scratch/bp72/rs9777/sim_results_roams/"

.libPaths("R_lib")
library(roams)
library(MASS)
library(stringr)
library(future)
library(furrr)
library(purrr)
library(dplyr)
library(tidyr)
library(tibble)

simulate_data_soft = function(
    sample_sizes = c(100, 200, 500, 1000),
    samples = 100,
    n_oos = 20,
    contamination = 0.1,
    distance = 5,
    trajectory_bias = 0.7,     # 0 = random direction, 1 = pure velocity direction
    multi_level_distances = c(distance - 2, distance, distance + 2),
    phi_coef = 0.8,
    sigma2_w_lon = 0.1,
    sigma2_w_lat = 0.1,
    sigma2_v_lon = 0.4,
    sigma2_v_lat = 0.4,
    initial_state = c(0, 0, 0, 0),
    seed = NA
) {

  if (length(n_oos) > 1 || length(contamination) > 1 ||
      length(distance) > 1 ||
      length(trajectory_bias) > 1 || trajectory_bias < 0 || trajectory_bias > 1 ||
      length(multi_level_distances) != 3 || length(phi_coef) > 1 ||
      length(sigma2_w_lon) > 1 || length(sigma2_w_lat) > 1 ||
      length(sigma2_v_lon) > 1 || length(sigma2_v_lat) > 1 ||
      length(initial_state) != 4) {
    stop("Check argument lengths. 'trajectory_bias' must be a single value in [0, 1], 'multi_level_distances' must be length 3, 'initial_state' must be length 4. All other arguments must be length 1.")
  }

  samples_seq       = 1:samples
  ns                = sample_sizes
  multi_level_dists = sort(multi_level_distances)
  dist              = distance
  cont              = contamination

  if (!is.na(seed)) set.seed(seed)

  data_sets = tibble::tibble(
    sample         = numeric(),
    n              = numeric(),
    setting        = character(),
    y_oos          = list(),
    x_oos          = list(),
    y              = list(),
    y_clean        = list(),
    x              = list(),
    outlier_locs   = list(),
    outlier_levels = list()
  )

  settings = c("clean data", "velocity-biased", "velocity multi-level")

  # DCRW matrices
  Phi = diag(c(1 + phi_coef, 1 + phi_coef, 0, 0))
  Phi[1, 3] = -phi_coef
  Phi[2, 4] = -phi_coef
  Phi[3, 1] = 1
  Phi[4, 2] = 1
  A  = diag(4)[1:2, ]
  Q  = diag(c(sigma2_w_lon, sigma2_w_lat, 0, 0))
  R  = diag(c(sigma2_v_lon, sigma2_v_lat))
  x0 = initial_state

  for (setting in settings) {
    for (n in ns) {
      for (sample_number in samples_seq) {

        # --- State process ---
        x = matrix(nrow = n + n_oos, ncol = 4)
        x[1, ] = Phi %*% x0 + MASS::mvrnorm(mu = rep(0, 4), Sigma = Q)
        for (t in 2:(n + n_oos)) {
          x[t, ] = Phi %*% x[t - 1, ] + MASS::mvrnorm(mu = rep(0, 4), Sigma = Q)
        }

        # --- Clean observations ---
        y_clean = matrix(nrow = n + n_oos, ncol = 2)
        for (t in 1:(n + n_oos)) {
          y_clean[t, ] = A %*% x[t, ] + MASS::mvrnorm(mu = rep(0, 2), Sigma = R)
        }

        y              = y_clean
        outlier_locs   = rep(0, n)
        outlier_levels = NA

        if (setting != "clean data") {

          # Outlier locations: first time point cannot be an outlier
          outlier_locs = c(
            0,
            sample(c(rep(1, round(n * cont)), rep(0, n - round(n * cont) - 1)))
          )

          # Velocity directions from the DCRW state:
          # x[t, 1:2] is current position, x[t, 3:4] is the lagged position (= x[t-1, 1:2]).
          # Their difference gives the true instantaneous velocity at each time point.
          vel_raw   = x[1:n, 1:2] - x[1:n, 3:4]         # n x 2
          vel_norms = sqrt(rowSums(vel_raw^2))
          vel_norms[vel_norms == 0] = 1                    # guard against zero velocity
          vel_dirs  = vel_raw / vel_norms                  # unit velocity directions

          # Random angles for the uninformative component
          rand_angles = runif(n, max = 2 * pi)
          vel_angles  = atan2(vel_dirs[, 2], vel_dirs[, 1])

          # Angular interpolation: rotate from the random angle toward the velocity angle
          # by a fraction trajectory_bias, taking the shortest arc (difference in [-pi, pi]).
          # At trajectory_bias = 0 the direction is purely random; at 1 it is purely velocity.
          angle_diff    = vel_angles - rand_angles
          angle_diff    = ((angle_diff + pi) %% (2 * pi)) - pi   # wrap to [-pi, pi]
          biased_angles = rand_angles + trajectory_bias * angle_diff
          biased_dirs   = cbind(cos(biased_angles), sin(biased_angles))  # n x 2 unit vectors

          if (setting == "velocity-biased") {
            # Outlier displaced at fixed 'distance' in the biased direction
            for (t in 1:n) {
              if (outlier_locs[t] == 1) {
                y[t, ] = y_clean[t, ] + dist * biased_dirs[t, ]
              }
            }

          } else if (setting == "velocity multi-level") {
            # Outlier displaced at a level-specific distance in the biased direction
            levels = sample(c(
              rep(multi_level_dists[1], floor(n * cont / 3)),
              rep(multi_level_dists[2], floor(n * cont / 3)),
              rep(multi_level_dists[3], n * cont - 2 * floor(n * cont / 3))
            ))

            outlier_levels = rep(0, n)
            counter = 1
            for (t in 1:n) {
              if (outlier_locs[t] == 1) {
                outlier_levels[t] = levels[counter]
                y[t, ] = y_clean[t, ] + outlier_levels[t] * biased_dirs[t, ]
                counter = counter + 1
              }
            }

          }
        }

        data_sets = data_sets |>
          tibble::add_row(
            sample         = sample_number,
            n              = n,
            setting        = setting,
            y_oos          = list(y[(n + 1):(n + n_oos), ]),
            x_oos          = list(x[(n + 1):(n + n_oos), ]),
            y              = list(y[1:n, ]),
            y_clean        = list(y_clean[1:n, ]),
            x              = list(x[1:n, ]),
            outlier_locs   = list(outlier_locs),
            outlier_levels = list(outlier_levels)
          )
      }
    }
  }

  data_sets |>
    dplyr::select(sample, n, setting, y, x, y_oos, x_oos, y_clean, outlier_locs, outlier_levels)
}

seed = 2025.4
true_par = c(0.8, 0.1, 0.1, 0.4, 0.4)
tbs = c(0, 0.25, 0.5, 0.75, 1)
data_study4 = tibble()

for (tb in tbs) {
  data_study4 = data_study4 |>
    bind_rows(
      simulate_data_soft(
        sample_sizes = c(100, 200, 500, 1000),
        samples = num_samples,
        n_oos = 20,
        contamination = 0.1,
        distance = 5,
        trajectory_bias = tb,
        phi_coef = true_par[1],
        sigma2_w_lon = true_par[2],
        sigma2_w_lat = true_par[3],
        sigma2_v_lon = true_par[4],
        sigma2_v_lat = true_par[5],
        initial_state = c(0, 0, 0, 0),
        seed = as.numeric(str_remove(seed, "seed"))
      ) |>
        mutate(trajectory_bias = tb)
    ) |>
    filter(setting == "velocity-biased")
}

# Parallelisation
plan(multisession, workers = num_cores)

# ROAMS --- SOFT
roams_start = Sys.time()
model_fits_roams_soft = data_study4 |>
  mutate(
    model_list = furrr::future_map(y, function(y) {

      x0 = c(y[1,], y[1,])

      build = function(par, x0) {

        phi_coef = par[1]
        Phi = diag(c(1+phi_coef, 1+phi_coef, 0, 0))
        Phi[1,3] = -phi_coef
        Phi[2,4] = -phi_coef
        Phi[3,1] = 1
        Phi[4,2] = 1

        A = diag(4)[1:2,]
        Q = diag(c(par[2], par[3], 0, 0))
        R = diag(c(par[4], par[5]))
        P0 = diag(rep(0, 4))

        specify_SSM(
          state_transition_matrix = Phi,
          state_noise_var = Q,
          observation_matrix = A,
          observation_noise_var = R,
          init_state_mean = x0,
          init_state_var = P0)
      }

      roams_SSM(
        y = y,
        init_par = c(0.5, rep(c(mad(diff(y[,1]), na.rm = TRUE)^2,
                                mad(diff(y[,2]), na.rm = TRUE)^2), 2)),
        build = build,
        build_args = list(x0 = x0),
        num_lambdas = 20,
        cores = 1,
        lower = c(0, rep(1e-12, 4)),
        upper = c(1, rep(Inf, 4)),
        B = 50,
        thresholding = "soft",
        lambda_min = 1
      )
    },
    .progress = TRUE),
    best_BIC_model = map(model_list, best_BIC_model),
    .keep = "none")
roams_end = Sys.time()
roams_end - roams_start
saveRDS(
  model_fits_roams_soft,
  file = paste0(data_path,
                "roams_soft_B",
                num_samples, "_",
                seed, ".rds")
)
rm(model_fits_roams_soft)  ## free up memory

# ROAMS --- HARD
roams_start = Sys.time()
model_fits_roams_hard = data_study4 |>
  mutate(
    model_list = furrr::future_map(y, function(y) {

      x0 = c(y[1,], y[1,])

      build = function(par, x0) {

        phi_coef = par[1]
        Phi = diag(c(1+phi_coef, 1+phi_coef, 0, 0))
        Phi[1,3] = -phi_coef
        Phi[2,4] = -phi_coef
        Phi[3,1] = 1
        Phi[4,2] = 1

        A = diag(4)[1:2,]
        Q = diag(c(par[2], par[3], 0, 0))
        R = diag(c(par[4], par[5]))
        P0 = diag(rep(0, 4))

        specify_SSM(
          state_transition_matrix = Phi,
          state_noise_var = Q,
          observation_matrix = A,
          observation_noise_var = R,
          init_state_mean = x0,
          init_state_var = P0)
      }

      roams_SSM(
        y = y,
        init_par = c(0.5, rep(c(mad(diff(y[,1]), na.rm = TRUE)^2,
                                mad(diff(y[,2]), na.rm = TRUE)^2), 2)),
        build = build,
        build_args = list(x0 = x0),
        num_lambdas = 20,
        cores = 1,
        lower = c(0, rep(1e-12, 4)),
        upper = c(1, rep(Inf, 4)),
        B = 50,
        thresholding = "hard",
        lambda_min = 2
      )
    },
    .progress = TRUE),
    best_BIC_model = map(model_list, best_BIC_model),
    .keep = "none")
roams_end = Sys.time()
roams_end - roams_start
saveRDS(
  model_fits_roams_hard,
  file = paste0(data_path,
                "roams_hard_B",
                num_samples, "_",
                seed, ".rds")
)
rm(model_fits_roams_hard)  ## free up memory

plan(sequential)




