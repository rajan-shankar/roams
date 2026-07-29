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


# Extract a sample and plot the 2D observations
sample_id = 3
y_mat = data_study1_n100 |>
  dplyr::filter(setting == "cluster", sample == sample_id) |>
  dplyr::pull(y) |>
  purrr::pluck(1)

ggplot2::ggplot(as.data.frame(y_mat), ggplot2::aes(x = V1, y = V2)) +
  ggplot2::geom_point()


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
  timepoint = 1:100,
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
b = attach_insample_info(res_ml$soft_list[[1]][[12]])
c = attach_insample_info(res_ml$soft_list[[1]][[15]])
d = attach_insample_info(res_ml$soft_list[[1]][[20]])

path_data = tibble(
  y1_BIC = (res_ml$y[[1]] - res_ml$soft_BIC[[1]]$gamma)[,1],
  y2_BIC = (res_ml$y[[1]] - res_ml$soft_BIC[[1]]$gamma)[,2],
  y1_10 = (res_ml$y[[1]] - res_ml$soft_list[[1]][[12]]$gamma)[,1],
  y2_10 = (res_ml$y[[1]] - res_ml$soft_list[[1]][[12]]$gamma)[,2],
  y1_20 = (res_ml$y[[1]] - res_ml$soft_list[[1]][[15]]$gamma)[,1],
  y2_20 = (res_ml$y[[1]] - res_ml$soft_list[[1]][[15]]$gamma)[,2],
  y1_30 = (res_ml$y[[1]] - res_ml$soft_list[[1]][[20]]$gamma)[,1],
  y2_30 = (res_ml$y[[1]] - res_ml$soft_list[[1]][[20]]$gamma)[,2],
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
  y1_30_fit = d$smoothed_observations[,1],
  y2_30_fit = d$smoothed_observations[,2],
  outlier_locs = factor(res_ml$outlier_locs[[1]])
)

p1 = path_data |>
  ggplot() +
  geom_segment(aes(x = y1, xend = y1_BIC, y = y2, yend = y2_BIC),
               colour = "orange", linetype = "dashed") +
  geom_point(aes(x = y1, y = y2, colour = outlier_locs),
             alpha = 0.4, size = 2) +
  geom_point(aes(x = y1_BIC, y = y2_BIC, colour = outlier_locs),
             size = 2) +
  geom_path(aes(x = x1, y = x2), colour = "black",
            linewidth = 0.8) +
  geom_path(aes(x = y1_BIC_fit, y = y2_BIC_fit), colour = "orange",
            linewidth = 0.8) +
  scale_color_manual(values = c("grey", "#D24644FF")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(size = 11, hjust = 0.5)) +
  labs(x = NULL, y = NULL, title = latex2exp::TeX("$\\lambda = 1.9$ (BIC-selected)"))

p2 = path_data |>
  ggplot() +
  geom_segment(aes(x = y1, xend = y1_10, y = y2, yend = y2_10),
               colour = "orange", linetype = "dashed") +
  geom_point(aes(x = y1, y = y2, colour = outlier_locs),
             alpha = 0.4, size = 2) +
  geom_point(aes(x = y1_10, y = y2_10, colour = outlier_locs),
             size = 2) +
  geom_path(aes(x = x1, y = x2), colour = "black",
            linewidth = 0.8) +
  geom_path(aes(x = y1_10_fit, y = y2_10_fit), colour = "orange",
            linewidth = 0.8) +
  scale_color_manual(values = c("grey", "#D24644FF")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(size = 11, hjust = 0.5)) +
  labs(x = NULL, y = NULL, title = latex2exp::TeX("$\\lambda = 3.4$"))

p3 = path_data |>
  ggplot() +
  geom_segment(aes(x = y1, xend = y1_20, y = y2, yend = y2_20),
               colour = "orange", linetype = "dashed") +
  geom_point(aes(x = y1, y = y2, colour = outlier_locs),
             alpha = 0.4, size = 2) +
  geom_point(aes(x = y1_20, y = y2_20, colour = outlier_locs)) +
  geom_path(aes(x = x1, y = x2), colour = "black",
            linewidth = 0.8) +
  geom_path(aes(x = y1_20_fit, y = y2_20_fit), colour = "orange",
            linewidth = 0.8) +
  scale_color_manual(values = c("grey", "#D24644FF")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(size = 11, hjust = 0.5)) +
  labs(x = NULL, y = NULL, title = latex2exp::TeX("$\\lambda = 4.1$"))

p4 = path_data |>
  ggplot() +
  geom_segment(aes(x = y1, xend = y1_30, y = y2, yend = y2_30),
               colour = "orange", linetype = "dashed") +
  geom_point(aes(x = y1, y = y2, colour = outlier_locs),
             alpha = 0.4, size = 2) +
  geom_point(aes(x = y1_30, y = y2_30, colour = outlier_locs),
             size = 2) +
  geom_path(aes(x = x1, y = x2), colour = "black",
            linewidth = 0.8) +
  geom_path(aes(x = y1_30_fit, y = y2_30_fit), colour = "orange",
            linewidth = 0.8) +
  scale_color_manual(values = c("grey", "#D24644FF")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none", plot.title = element_text(size = 11, hjust = 0.5)) +
  labs(x = NULL, y = NULL, title = latex2exp::TeX("$\\lambda = 5.2$"))

# Shared legend
leg_features = c("clean observation", "true outlier", "true path", "ROAMS path", "mean-shift")
leg_colours  = c("grey",             "#D24644FF",    "black",      "orange",     "orange")
leg_labels   = c("Clean observation","True outlier", "True path",  "ROAMS path", "Estimated mean-shift")

ld = data.frame(x = c(0, 1), y = c(0, 1),
                feature = rep(leg_features, each = 2))

legend_plot = ggplot(ld, aes(x = x, y = y, colour = feature)) +
  geom_point() +
  geom_line() +
  scale_colour_manual(
    values   = setNames(leg_colours,  leg_features),
    labels   = setNames(leg_labels,   leg_features),
    breaks   = leg_features
  ) +
  theme_void(base_size = 14) +
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.key.width = unit(1.2, "cm")) +
  guides(colour = guide_legend(
    nrow = 2, byrow = TRUE,
    override.aes = list(
      shape     = c(16,      16,           NA,       NA,        NA),
      linetype  = c("blank", "blank",      "solid",  "solid",   "dashed"),
      alpha     = c(1,       1,            1,        1,         1),
      linewidth = c(0,       0,            0.8,      0.8,       0.8)
    )
  ))

legend_obj = cowplot::get_legend(legend_plot)

(p4 + p3 + p2 + p1 + plot_layout(nrow = 1, axes = "collect")) /
  patchwork::wrap_elements(legend_obj) +
  plot_layout(heights = c(4, 1))

ggsave(paste0("soft_thresholding/soft_thresholding_cluster.pdf"),
       width = 7, height = 3)


# Mean-shift vector norm path plot
mean_shifts = get_attribute(res_ml$soft_list[[1]], "gamma")
lambda_vals = get_attribute(res_ml$soft_list[[1]], "lambda")
true_outliers = res_ml$outlier_locs[[1]]

path_df = purrr::map2_dfr(mean_shifts, lambda_vals, function(gamma, lam) {
  tibble::tibble(
    timepoint    = 1:nrow(gamma),
    lambda       = lam,
    norm         = sqrt(gamma[, 1]^2 + gamma[, 2]^2),
    true_outlier = as.logical(true_outliers)
  )
})

# Restrict to timepoints with at least one non-zero norm across the lambda path
active_timepoints = path_df |>
  dplyr::group_by(timepoint) |>
  dplyr::summarise(ever_nonzero = any(norm > 0)) |>
  dplyr::filter(ever_nonzero) |>
  dplyr::pull(timepoint)

# Hard thresholding path data
mean_shifts_hard = get_attribute(res_ml$hard_list[[1]], "gamma")
lambda_vals_hard = get_attribute(res_ml$hard_list[[1]], "lambda")

path_df_hard = purrr::map2_dfr(mean_shifts_hard, lambda_vals_hard, function(gamma, lam) {
  tibble::tibble(
    timepoint    = 1:nrow(gamma),
    lambda       = lam,
    norm         = sqrt(gamma[, 1]^2 + gamma[, 2]^2),
    true_outlier = as.logical(true_outliers)
  )
})

active_timepoints_hard = path_df_hard |>
  dplyr::group_by(timepoint) |>
  dplyr::summarise(ever_nonzero = any(norm > 0)) |>
  dplyr::filter(ever_nonzero) |>
  dplyr::pull(timepoint)

# Shared y-axis limit for fair comparison
y_max = max(
  path_df      |> dplyr::filter(timepoint %in% active_timepoints)      |> dplyr::pull(norm),
  path_df_hard |> dplyr::filter(timepoint %in% active_timepoints_hard) |> dplyr::pull(norm)
)

make_norm_plot = function(df, active_tp, title_str, bic_lambda) {
  df |>
    dplyr::filter(timepoint %in% active_tp) |>
    ggplot(aes(x = lambda, y = norm, group = timepoint, colour = true_outlier)) +
    geom_vline(xintercept = bic_lambda, linetype = "dashed", colour = "black") +
    geom_line(alpha = 0.7) +
    scale_x_reverse() +
    scale_y_continuous(limits = c(0, y_max)) +
    scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "#D24644FF"),
                        labels = c("FALSE" = "Clean observation", "TRUE" = "True outlier")) +
    labs(x = latex2exp::TeX("$\\lambda$"),
         y = latex2exp::TeX("Mean-shift vector norm"),
         colour = NULL,
         subtitle = title_str) +
    theme_bw(base_size = 14) +
    theme(panel.grid                = element_blank(),
          legend.position           = "bottom",
          plot.subtitle             = element_text(hjust = 0.5))
}

p_soft = make_norm_plot(path_df,      active_timepoints,      "Soft thresholding",
                        bic_lambda = res_ml$soft_BIC[[1]]$lambda)
p_hard = make_norm_plot(path_df_hard, active_timepoints_hard, "Hard thresholding",
                        bic_lambda = res_ml$hard_BIC[[1]]$lambda)

p_soft + p_hard +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.position = "bottom")

ggsave("soft_thresholding/soft_hard_mean_shift.pdf",
       width = 7, height = 3.5)

# Hard vs soft gamma scatter plot
roams_soft_B5_2025.4 <- readRDS("soft_thresholding/roams_soft_B5_2025.4.rds")
roams_hard_B5_2025.4 <- readRDS("soft_thresholding/roams_hard_B5_2025.4.rds")
colnames(roams_soft_B5_2025.4) = c("soft_list", "soft_BIC")
colnames(roams_hard_B5_2025.4) = c("hard_list", "hard_BIC")
res = bind_cols(data_study4, roams_hard_B5_2025.4, roams_soft_B5_2025.4)
example_vb = res |> filter(sample == 5, n == 1000, trajectory_bias == 1)
example_notvb = res |> filter(sample == 5, n == 1000, trajectory_bias == 0)

scatter_df = bind_rows(
  tibble(
    hard_norm    = sqrt(rowSums(example_notvb$hard_BIC[[1]]$gamma^2)),
    soft_norm    = sqrt(rowSums(example_notvb$soft_BIC[[1]]$gamma^2)),
    true_outlier = as.logical(example_notvb$outlier_locs[[1]]),
    setting      = 'paste("Random-direction outliers  (", eta==0, ")")'
  ),
  tibble(
    hard_norm    = sqrt(rowSums(example_vb$hard_BIC[[1]]$gamma^2)),
    soft_norm    = sqrt(rowSums(example_vb$soft_BIC[[1]]$gamma^2)),
    true_outlier = as.logical(example_vb$outlier_locs[[1]]),
    setting      = 'paste("Velocity-biased outliers  (", eta==1, ")")'
  )
)

ggplot(scatter_df, aes(x = soft_norm, y = hard_norm, colour = true_outlier)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_point(alpha = 0.7, size = 2) +
  scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "#D24644FF"),
                      labels = c("FALSE" = "Clean observation", "TRUE" = "True outlier")) +
  facet_wrap(~ setting, scales = "fixed", labeller = label_parsed) +
  labs(x = "Mean-shift norm (soft)",
       y = "Mean-shift norm (hard)",
       colour = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(paste0("soft_thresholding/hard_soft_shift.pdf"),
       width = 6, height = 3)

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

