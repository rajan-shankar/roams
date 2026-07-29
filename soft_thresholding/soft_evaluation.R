# Evals
roams_hard = read_rds("soft_thresholding/roams_hard_B100_2025.4.rds") |> dplyr::select(best_BIC_model)
roams_soft = read_rds("soft_thresholding/roams_soft_B100_2025.4.rds") |> dplyr::select(best_BIC_model)
model_fits = bind_cols(
  data_study4,
  roams_hard |> select(best_BIC_model) |> rename(roams_hard = best_BIC_model),
  roams_soft |> select(best_BIC_model) |> rename(roams_soft = best_BIC_model)
)

source("vignettes/scripts/evaluation_functions.R")
get_MSFE = function(model, method, y_oos) {

  x0_oos = c(y_oos[1,], y_oos[1,])
  build_oos = function(par, x0_oos) {

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
      init_state_mean = x0_oos,
      init_state_var = P0)
  }

  output = oos_filter(y_oos, model, build_oos,
                      build_args = list(x0_oos = x0_oos),
                      multiplier = 1, threshold = Inf)

  pred = output$predicted_observations
  SFEs = rowSums((pred - y_oos)^2)
  MSFE = mean(SFEs)
  return(MSFE)
}
get_insample_fit_error = function(model, x) {
  model_extra = attach_insample_info(model)
  fit = model_extra$predicted_observations
  squared_errors = rowSums((fit - x[,1:2])^2)
  RMSE = sqrt(mean(squared_errors))
  return(RMSE)
}


row_num = 400
rs = roams_soft |>
  slice(row_num) |>
  pull(best_BIC_model) |>
  pluck(1)

rh = roams_hard |>
  slice(row_num) |>
  pull(best_BIC_model) |>
  pluck(1)

# Plot

tibble(soft_gamma = sqrt(rowSums((rs$gamma)^2)),
       hard_gamma = sqrt(rowSums((rh$gamma)^2)),
       true_outliers = data_study4 |> slice(row_num) |> pull(outlier_locs) |> pluck(1)) |>
ggplot() +
  aes(x = soft_gamma, y = hard_gamma, colour = factor(true_outliers)) +
  geom_point() +
  geom_abline()

res_ml = data_study4 |> slice(row_num)
path_data = tibble(
  y1 = res_ml$y[[1]][,1],
  y2 = res_ml$y[[1]][,2],
  x1 = res_ml$x[[1]][,1],
  x2 = res_ml$x[[1]][,2],
  outlier_locs = factor(res_ml$outlier_locs[[1]]),
  y1_hard_fit = attach_insample_info(rh)$smoothed_observations[,1],
  y2_hard_fit = attach_insample_info(rh)$smoothed_observations[,2],
  flagged = (sqrt(rowSums((rh$gamma)^2)) != 0)
)

path_data |>
  ggplot() +
  geom_point(aes(x = y1, y = y2, colour = outlier_locs)) +
  geom_point(data = path_data |>
               filter(flagged),
             aes(x = y1, y = y2),
             pch = 1, colour = "orange", size = 3) +
  geom_path(aes(x = x1, y = x2), colour = "black",
            linewidth = 0.8) +
  geom_path(aes(x = y1_hard_fit, y = y2_hard_fit), colour = "orange",
            linewidth = 0.8) +
  scale_color_manual(values = c("darkgrey", "red")) +
  theme_bw() +
  theme(legend.position = "none")
#



eval_start = Sys.time()
eval_metrics = model_fits |>
  pivot_longer(starts_with("roams_"),
               names_to = "method",
               values_to = "model") |>
  mutate(
    MSFE = pmap_dbl(list(model, method, y_oos), .f = get_MSFE),
    fit_error = pmap_dbl(list(model, x), .f = get_insample_fit_error),
    TPR = map2_dbl(model, outlier_locs, .f = get_TPR),
    TNR = map2_dbl(model, outlier_locs, .f = get_TNR),
    par_dist_phi       = map_dbl(
      model, .f = ~ get_par_dist_phi(.x, true_par)),
    par_dist_obs_var   = map_dbl(
      model, .f = ~ get_par_dist_obs_var(.x, true_par)),
    par_dist_state_var = map_dbl(
      model, .f = ~ get_par_dist_state_var(.x, true_par)),
    iterations = map2_dbl(model, method, .f = ~ ifelse(
      str_starts(.y, "roams"),
      .x$iterations,
      NA
    )),
    lambda = map2_dbl(model, method, .f = ~ ifelse(
      str_starts(.y, "roams"),
      .x$lambda,
      NA
    ))
  ) |>
  select(-c(model,
            y_oos, x_oos, y, y_clean, x,
            outlier_locs, outlier_levels))  ## reduce size of object by removing data set information

eval_end = Sys.time()
eval_end - eval_start


theme_set(theme_bw(base_size = 14))
eval_metrics |>
  group_by(n, method, trajectory_bias) |>
  summarise(mean = mean(MSFE),
            se = sd(MSFE) / sqrt(n())) |>
  ggplot() +
  aes(x = trajectory_bias, y = mean, colour = method,
      group = method, ymin = mean - se, ymax = mean + se) +
  geom_point(size = 2) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(fill = method), alpha = 0.2, colour = NA) +
  facet_wrap(~ factor(n), nrow = 1,
             labeller = as_labeller(function(x) paste0("n = ", x))) +
  scale_colour_discrete(labels = c("roams_hard" = "ROAMS (hard thresholding)",
                                   "roams_soft" = "ROAMS (soft thresholding)")) +
  scale_fill_discrete(labels = c("roams_hard" = "ROAMS (hard thresholding)",
                                 "roams_soft" = "ROAMS (soft thresholding)")) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", "0.25", "0.5", "0.75", "1")) +
  # scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "bottom") +
  labs(y = "MSFE",
       x = "Velocity bias", colour = NULL, fill = NULL)

ggsave(paste0("soft_thresholding/soft_MSFE.pdf"),
       width = 7, height = 3)
# ggsave(paste0("soft_thresholding/soft_insample_forecast_error.pdf"),
#        width = 7, height = 3)

