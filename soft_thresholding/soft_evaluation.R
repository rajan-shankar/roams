# Evals
roams_hard = read_rds("soft_thresholding/roams_hard_B5_2025.4.rds")
roams_soft = read_rds("soft_thresholding/roams_soft_B5_2025.4.rds")
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
  fit = model_extra$smoothed_observations
  squared_errors = rowSums((fit - x[,1:2])^2)
  RMSE = sqrt(mean(squared_errors))
  return(RMSE)
}

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

eval_metrics |>
  group_by(n, method, trajectory_bias) |>
  summarise(metric = mean(fit_error)) |>
  ggplot() +
  aes(x = trajectory_bias, y = metric, colour = method,
      group = method) +
  geom_point() +
  geom_line() +
  facet_wrap(~ factor(n))


