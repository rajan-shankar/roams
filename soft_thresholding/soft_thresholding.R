library(tidyverse)
library(patchwork)

roams_start = Sys.time()
model_fits_roams_study1 = data_study1 |>
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

        # x0 = c(y[1,], y[1,])
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

res_ml = res |>
  filter(setting == "cluster")

res_ml$hard_list[[1]] |> autoplot()
res_ml$hard_list[[1]] |> autoplot(attribute1 = "iterations")
res_ml$hard_BIC[[1]]

res_ml$soft_list[[1]] |> autoplot()
res_ml$soft_list[[1]] |> autoplot(attribute1 = "iterations")
res_ml$soft_BIC[[1]]

gamma_data = tibble(
  hard = apply(res_ml$hard_BIC[[1]]$gamma, 1, function(x) sqrt(sum(x^2))),
  soft = apply(res_ml$soft_BIC[[1]]$gamma, 1, function(x) sqrt(sum(x^2))),
  timepoint = 1:200,
  outlier_levels = res_ml$outlier_levels[[1]],
  outlier_locs = factor(res_ml$outlier_locs[[1]])
)

gamma_data |>
  ggplot(aes(x = hard, y = soft, colour = outlier_locs)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

gamma_data |>
  ggplot(aes(x = timepoint)) +
  geom_segment(aes(y = 0, yend = hard),
               linewidth = 1) +
  geom_segment(aes(y = 0, yend = soft),
               color = "red",
               linewidth = 1) +
  geom_point(aes(y = outlier_levels),
             color = "blue")

a = attach_insample_info(res_ml$soft_BIC[[1]])
b = attach_insample_info(res_ml$soft_list[[1]][[10]])
c = attach_insample_info(res_ml$soft_list[[1]][[20]])

path_data = tibble(
  y1_BIC = (res_ml$y[[1]] - res_ml$soft_BIC[[1]]$gamma)[,1],
  y2_BIC = (res_ml$y[[1]] - res_ml$soft_BIC[[1]]$gamma)[,2],
  y1_10 = (res_ml$y[[1]] - res_ml$soft_list[[1]][[10]]$gamma)[,1],
  y2_10 = (res_ml$y[[1]] - res_ml$soft_list[[1]][[10]]$gamma)[,2],
  y1_20 = (res_ml$y[[1]] - res_ml$soft_list[[1]][[20]]$gamma)[,1],
  y2_20 = (res_ml$y[[1]] - res_ml$soft_list[[1]][[20]]$gamma)[,2],
  y1 = res_ml$y[[1]][,1],
  y2 = res_ml$y[[1]][,2],
  x1 = res_ml$x[[1]][,1],
  x2 = res_ml$x[[1]][,2],
  y1_BIC_fit = a$smoothed_observations[,1],
  y2_BIC_fit = a$smoothed_observations[,2],
  y1_10_fit = b$smoothed_observations[,1],
  y2_10_fit = b$smoothed_observations[,2],
  y1_20_fit = c$smoothed_observations[,1],
  y2_20_fit = c$smoothed_observations[,2],
  outlier_locs = factor(res_ml$outlier_locs[[1]])
)

p1 = path_data |>
  ggplot() +
  geom_segment(aes(x = y1, xend = y1_BIC, y = y2, yend = y2_BIC),
               colour = "orange", linetype = "dashed") +
  geom_point(aes(x = y1, y = y2, colour = outlier_locs), alpha = 0.4) +
  geom_point(aes(x = y1_BIC, y = y2_BIC, colour = outlier_locs)) +
  geom_path(aes(x = x1, y = x2), colour = "black",
            linewidth = 0.8) +
  geom_path(aes(x = y1_BIC_fit, y = y2_BIC_fit), colour = "orange",
            linewidth = 0.8) +
  scale_color_manual(values = c("darkgrey", "red")) +
  theme_bw() +
  theme(legend.position = "none")

p2 = path_data |>
  ggplot() +
  geom_segment(aes(x = y1, xend = y1_10, y = y2, yend = y2_10),
               colour = "orange", linetype = "dashed") +
  geom_point(aes(x = y1, y = y2, colour = outlier_locs), alpha = 0.4) +
  geom_point(aes(x = y1_10, y = y2_10, colour = outlier_locs)) +
  geom_path(aes(x = x1, y = x2), colour = "black",
            linewidth = 0.8) +
  geom_path(aes(x = y1_10_fit, y = y2_10_fit), colour = "orange",
            linewidth = 0.8) +
  scale_color_manual(values = c("darkgrey", "red")) +
  theme_bw() +
  theme(legend.position = "none")

p3 = path_data |>
  ggplot() +
  geom_segment(aes(x = y1, xend = y1_20, y = y2, yend = y2_20),
               colour = "orange", linetype = "dashed") +
  geom_point(aes(x = y1, y = y2, colour = outlier_locs), alpha = 0.4) +
  geom_point(aes(x = y1_20, y = y2_20, colour = outlier_locs)) +
  geom_path(aes(x = x1, y = x2), colour = "black",
            linewidth = 0.8) +
  geom_path(aes(x = y1_20_fit, y = y2_20_fit), colour = "orange",
            linewidth = 0.8) +
  scale_color_manual(values = c("darkgrey", "red")) +
  theme_bw() +
  theme(legend.position = "none")

p3 + p2 + p1 + plot_layout(guides = "collect")



# sim_data = simulate_data_soft(
#   sample_sizes = 100,
#   samples = 1,
#   trajectory_bias = 1,
#   seed = 1
# ) |>
#   filter(setting == "velocity-biased")
#
# sim_data |>
#   select(sample:x, outlier_locs) |>
#   mutate(y1 = map(y, ~ .[,1]),
#          y2 = map(y, ~ .[,2]),
#          x1 = map(x, ~ .[,1]),
#          x2 = map(x, ~ .[,2])) |>
#   select(-y, -x) |>
#   unnest(everything()) |>
#   mutate(outlier_locs = factor(outlier_locs)) |>
#   ggplot() +
#   geom_segment(data = . %>% filter(outlier_locs == 1),
#                aes(x = x1, xend = y1, y = x2, yend = y2),
#                colour = "red", linetype = "dashed") +
#   geom_point(aes(x = y1, y = y2, colour = outlier_locs)) +
#   geom_path(aes(x = x1, y = x2)) +
#   facet_wrap(~ setting + sample, nrow = 3, scales= "free") +
#   scale_color_manual(values = c("darkgrey", "red")) +
#   theme_bw() +
#   theme(legend.position = "none")

