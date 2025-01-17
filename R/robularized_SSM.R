#' A Cat Function
#'
#' This function allows you to express your love of cats.
#'
#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @param love Do you love cats? Defaults to TRUE.
#' @details
#' Additional details...
#' @returns description
#' @export
robularized_SSM = function(
    y,
    init_par,
    build,
    num_lambdas = 10,
    lowest_lambda = 2,
    highest_lambda = NA,
    num_lambdas_crowding = 0,
    cores = 1,
    B = 20,
    lower = NA,
    upper = NA
    ) {

  # Classical fit by using large lambda
  classical = run_IPOD(y = y,
                       lambda = 100,
                       init_par = init_par,
                       build = build,
                       B = B,
                       lower = lower,
                       upper = upper)

  # Highest lambda is the supremum norm of mahalanobis residuals of classical fit
  if (is.na(highest_lambda)) {
    highest_lambda = max(fn_filter(classical$par,
                                   y,
                                   classical$gamma,
                                   build)$mahalanobis_residuals
                         )
  }

  # lambda grid
  lambdas = seq(lowest_lambda,
                highest_lambda,
                length.out = num_lambdas)

  # Fit models across the grid
  model_list = lambda_grid(y = y,
                           lambdas = lambdas,
                           #init_par = classical$par,
                           init_par = init_par,
                           build = build,
                           cores = cores,
                           lower = lower,
                           upper = upper,
                           B = B)

  if (num_lambdas_crowding != 0) {
    BIC = get_attribute(model_list, "BIC")
    prop_outlying = get_attribute(model_list, "prop_outlying")

    # Take BIC diffs and set any that are left of prop_outlying >= 0.5 to 0
    bad_up_to_and_including = sum(prop_outlying >= 0.5)
    diffs = abs(diff(BIC))
    if (bad_up_to_and_including >= 2) {
      diffs[1:(bad_up_to_and_including - 1)] = 0
    }

    # Find largest BIC diff
    diff_argmax = which.max(diffs)

    # Crowding lambda grid
    crowding_lambdas = seq(lambdas[diff_argmax],
                           lambdas[diff_argmax + 1],
                           length.out = num_lambdas_crowding + 2)[2:(num_lambdas_crowding+1)]

    # Fit models across the crowding grid
    model_crowding_list = lambda_grid(y = y,
                                      lambdas = crowding_lambdas,
                                      #init_par = classical$par,
                                      init_par = init_par,
                                      build = build,
                                      cores = cores,
                                      lower = lower,
                                      upper = upper,
                                      B = B)

    # Append to model_crowding_list to model_list at appropriate location
    model_list = append(model_list, model_crowding_list, after = diff_argmax)

    # Attach class
    class(model_list) = "robularized_SSM_list"
  }

  return(model_list)
}

lambda_grid = function(
    y,
    lambdas,
    init_par,
    build,
    cores,
    lower,
    upper,
    B
    ) {

  cl = parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl,
                                 export = list(fn_filter = fn_filter,
                                               run_IPOD = run_IPOD))
  model_list = foreach::foreach(
    i = 1:length(lambdas)) %dopar% {
                  IPOD_output = run_IPOD(y,
                                         lambdas[i],
                                         init_par,
                                         build,
                                         B,
                                         lower,
                                         upper)
                  filter_output = fn_filter(IPOD_output$par,
                                            y,
                                            IPOD_output$gamma,
                                            build)

                  model = c(IPOD_output, filter_output)
                  class(model) = "robularized_SSM"
                  return(model)
                  }
  parallel::stopCluster(cl)

  class(model_list) = "robularized_SSM_list"
  return(model_list)
}

run_IPOD = function(
    y,
    lambda,
    init_par,
    build,
    B,
    lower,
    upper
    ) {

  if (is.na(lower)[1]) {lower = rep(-Inf, length(init_par))}
  if (is.na(upper)[1]) {upper = rep(Inf, length(init_par))}

  n = ncol(y)
  dim_obs = nrow(y)
  par = init_par
  gamma = matrix(0, nrow = dim_obs, ncol = n)
  r = NA
  theta_old = par

  for (j in 1:B) {
    if (j != 1) {theta_old = res$par}
    res = stats::optim(
      par = par,
      fn = fn_filter,
      y = y,
      gamma = gamma,
      build = build,
      return_obj = TRUE,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper
      )

    if ((sum(res$par == lower) + sum(res$par == upper)) == 0) {
      par = res$par
    } else {
      par = init_par
    }

    filter_output = fn_filter(res$par, y, gamma, build)
    r = y - filter_output$predicted_observations
    gamma_old = gamma
    gamma = matrix(0, nrow = dim_obs, ncol = n)
    gamma[,filter_output$mahalanobis_residuals > lambda] = r[,filter_output$mahalanobis_residuals > lambda]
    gap = max(abs(gamma - gamma_old))
    gap_theta = max(abs(res$par - theta_old))

    nz = sum(colSums(abs(gamma_old)) != 0)
    prop_outlying = nz / n

    if ((gap < 1e-4) && (gap_theta < 1e-4)) {
      break
    }
    # new termination criterion
    if (prop_outlying >= 0.5) {
      break
    }
  }

  p = length(init_par)
  RSS = sum((r - gamma_old)^2, na.rm = TRUE)
  BIC = (n-p)*log(RSS/(n-p)) + (nz+1)*(log(n-p) + 1)
  #negloglik = n*fn_filter(model$par, gamma = gamma_old, y = y, return_obj = TRUE)

  return(list(
    "lambda" = lambda,
    "par" = res$par,
    "prop_outlying" = prop_outlying,
    "BIC" = BIC,
    "RSS" = RSS,
    "gamma" = gamma_old,
    "iterations" = j
    ))
}


# Missing value handling + smoother
fn_filter = function(
    par,
    y,
    gamma,
    build,
    return_obj = FALSE
) {

  SSM_specs = build(par)

  Phi = SSM_specs$state_transition_matrix
  Sigma_w = SSM_specs$state_noise_var
  A = SSM_specs$observation_matrix
  Sigma_v = SSM_specs$observation_noise_var
  x_tt = SSM_specs$init_state_mean
  P_tt = SSM_specs$init_state_var

  n = ncol(y)
  dim_obs = nrow(y)
  dim_state = nrow(Phi)

  x_tt_1 = NA
  P_tt_1 = NA
  y_tt_1 = NA
  S_t = NA
  objective = 0

  filtered_states = matrix(0, nrow = dim_state, ncol = n)
  filtered_observations = matrix(0, nrow = dim_obs, ncol = n)
  predicted_states = matrix(0, nrow = dim_state, ncol = n)
  predicted_observations = matrix(0, nrow = dim_obs, ncol = n)
  predicted_observations_var = list()
  filtered_states_var = list()  # V
  predicted_states_var = list()  # P
  mahalanobis_residuals = NA

  for (t in 1:n) {
    x_tt_1 = Phi %*% x_tt
    P_tt_1 = Phi %*% P_tt %*% t(Phi) + Sigma_w
    y_tt_1 = A %*% x_tt_1
    S_t = A %*% P_tt_1 %*% t(A) + Sigma_v
    inv_S_t = solve(S_t)

    if (is.na(y[1,t])) {
      x_tt = x_tt_1
      P_tt = P_tt_1

      filtered_states[,t] = x_tt
      filtered_observations[,t] = A %*% x_tt
      predicted_states[,t] = x_tt_1
      predicted_observations[,t] = y_tt_1
      predicted_observations_var[[t]] = S_t
      filtered_states_var[[t]] = P_tt
      predicted_states_var[[t]] = P_tt_1

    } else {

      if (sum(abs(gamma[,t])) == 0) {
        K_t = P_tt_1 %*% t(A) %*% inv_S_t
        x_tt = x_tt_1 + K_t %*% (y[,t] - y_tt_1)
        P_tt = P_tt_1 - K_t %*% A %*% P_tt_1
      } else {
        x_tt = x_tt_1
        P_tt = P_tt_1
      }

      filtered_states[,t] = x_tt
      filtered_observations[,t] = A %*% x_tt
      predicted_states[,t] = x_tt_1
      predicted_observations[,t] = y_tt_1
      predicted_observations_var[[t]] = S_t
      filtered_states_var[[t]] = P_tt
      predicted_states_var[[t]] = P_tt_1
    }
  }

  smoothed_states = matrix(0, nrow = dim_state, ncol = n)
  smoothed_states[,n] = filtered_states[,n]
  smoothed_states_var = list()
  smoothed_states_var[[n]] = filtered_states_var[[n]]
  smoothed_observations = matrix(0, nrow = dim_obs, ncol = n)
  smoothed_observations[,n] = A %*% smoothed_states[,n]
  for (t in (n-1):1) {

    C_t = filtered_states_var[[t]] %*% t(Phi) %*% solve(predicted_states_var[[t+1]])

    # if ((sum(abs(gamma[,t])) != 0) | is.na(y[1,t])) {
    #   smoothed_states[,t] = smoothed_states[,t+1]
    #   smoothed_states_var[[t]] = smoothed_states_var[[t+1]]
    # } else {
    smoothed_states[,t] = filtered_states[,t] + C_t %*% (smoothed_states[,t+1] - Phi %*% filtered_states[,t])
    smoothed_states_var[[t]] = filtered_states_var[[t]] + C_t %*% (smoothed_states_var[[t+1]] - predicted_states_var[[t+1]]) %*% t(C_t)
    # }

    smoothed_observations[,t] = A %*% smoothed_states[,t]
  }

  predicted_observations_updated = matrix(0, nrow = dim_obs, ncol = n)
  predicted_observations_var_updated = list()
  x_tt = SSM_specs$init_state_mean
  P_tt = SSM_specs$init_state_var
  for (t in 1:n) {
    x_tt_1 = Phi %*% x_tt
    P_tt_1 = Phi %*% P_tt %*% t(Phi) + Sigma_w
    y_tt_1 = A %*% x_tt_1
    S_t = A %*% P_tt_1 %*% t(A) + Sigma_v
    inv_S_t = solve(S_t)

    x_tt = smoothed_states[,t]
    P_tt = smoothed_states_var[[t]]

    predicted_observations_updated[,t] = y_tt_1
    predicted_observations_var_updated[[t]] = S_t

    if (!is.na(y[1,t])) {
      objective = objective + 1/(2*n) * ((sum(abs(gamma[,t])) == 0) * log(det(S_t)) + t(y[,t] - y_tt_1 - gamma[,t]) %*% inv_S_t %*% (y[,t] - y_tt_1 - gamma[,t]))
      mahalanobis_residuals[t] = drop(sqrt(t(y[,t] - y_tt_1) %*% inv_S_t %*% (y[,t] - y_tt_1)))
    } else {
      mahalanobis_residuals[t] = 0
    }
  }

  squared_errors = (y - predicted_observations_updated)^2

  if (return_obj) {
    return(objective)
  } else {
    return(list(
      "filtered_states" = filtered_states,
      "filtered_observations" = filtered_observations,
      "smoothed_states" = smoothed_states,
      "smoothed_observations" = smoothed_observations,
      "predicted_states" = predicted_states,
      "predicted_observations" = predicted_observations_updated,
      "predicted_observations_var" = predicted_observations_var_updated,
      "mahalanobis_residuals" = mahalanobis_residuals,
      "squared_errors" = squared_errors,

      "smoothed_states_var" = smoothed_states_var,
      "predicted_states_var" = predicted_states_var,
      "filtered_states_var" = filtered_states_var
    ))
  }
}



# Missing value handling
# fn_filter = function(
#     par,
#     y,
#     gamma,
#     build,
#     return_obj = FALSE
# ) {
#
#   SSM_specs = build(par)
#
#   Phi = SSM_specs$state_transition_matrix
#   Sigma_w = SSM_specs$state_noise_var
#   A = SSM_specs$observation_matrix
#   Sigma_v = SSM_specs$observation_noise_var
#   x_tt = SSM_specs$init_state_mean
#   P_tt = SSM_specs$init_state_var
#
#   n = ncol(y)
#   dim_obs = nrow(y)
#   dim_state = nrow(Phi)
#
#   x_tt_1 = NA
#   P_tt_1 = NA
#   y_tt_1 = NA
#   S_t = NA
#   objective = 0
#
#   if (!return_obj) {
#     filtered_states = matrix(0, nrow = dim_state, ncol = n)
#     filtered_observations = matrix(0, nrow = dim_obs, ncol = n)
#     predicted_states = matrix(0, nrow = dim_state, ncol = n)
#     predicted_observations = matrix(0, nrow = dim_obs, ncol = n)
#     predicted_observations_var = list()
#     mahalanobis_residuals = NA
#   }
#
#   for (t in 1:n) {
#     x_tt_1 = Phi %*% x_tt
#     P_tt_1 = Phi %*% P_tt %*% t(Phi) + Sigma_w
#     y_tt_1 = A %*% x_tt_1
#     S_t = A %*% P_tt_1 %*% t(A) + Sigma_v
#     inv_S_t = solve(S_t)
#
#     if (is.na(y[1,t])) {
#       x_tt = x_tt_1
#       P_tt = P_tt_1
#
#       if (return_obj) {
#         objective = objective + 0
#       } else {
#         filtered_states[,t] = x_tt
#         filtered_observations[,t] = A %*% x_tt
#         predicted_states[,t] = x_tt_1
#         predicted_observations[,t] = y_tt_1
#         predicted_observations_var[[t]] = S_t
#         mahalanobis_residuals[t] = 0
#       }
#
#     } else {
#
#       if (sum(abs(gamma[,t])) == 0) {
#         K_t = P_tt_1 %*% t(A) %*% inv_S_t
#         x_tt = x_tt_1 + K_t %*% (y[,t] - y_tt_1)
#         P_tt = P_tt_1 - K_t %*% A %*% P_tt_1
#       } else {
#         x_tt = x_tt_1
#         P_tt = P_tt_1
#       }
#       if (return_obj) {
#         objective = objective + 1/(2*n) * ((sum(abs(gamma[,t])) == 0) * log(det(S_t)) + t(y[,t] - y_tt_1 - gamma[,t]) %*% inv_S_t %*% (y[,t] - y_tt_1 - gamma[,t]))
#       } else {
#         filtered_states[,t] = x_tt
#         filtered_observations[,t] = A %*% x_tt
#         predicted_states[,t] = x_tt_1
#         predicted_observations[,t] = y_tt_1
#         predicted_observations_var[[t]] = S_t
#         mahalanobis_residuals[t] = drop(sqrt(t(y[,t] - y_tt_1) %*% inv_S_t %*% (y[,t] - y_tt_1)))
#       }
#     }
#
#   }
#
#   if (return_obj) {
#     return(objective)
#   } else {
#     return(list(
#       "filtered_states" = filtered_states,
#       "filtered_observations" = filtered_observations,
#       "predicted_states" = predicted_states,
#       "predicted_observations" = predicted_observations,
#       "predicted_observations_var" = predicted_observations_var,
#       "mahalanobis_residuals" = mahalanobis_residuals
#     ))
#   }
# }


# Handles exogenous inputs in A matrix
# fn_filter = function(
#     par,
#     y,
#     gamma,
#     build,
#     return_obj = FALSE
# ) {
#
#   SSM_specs = build(par)
#
#   Phi = SSM_specs$state_transition_matrix
#   Sigma_w = SSM_specs$state_noise_var
#   A = SSM_specs$observation_matrix
#   Sigma_v = SSM_specs$observation_noise_var
#   x_tt = SSM_specs$init_state_mean
#   P_tt = SSM_specs$init_state_var
#
#   n = ncol(y)
#   dim_obs = nrow(y)
#   dim_state = nrow(Phi)
#
#   if (class(A) != "list") {
#     A_mat = A
#     A = list()
#     for (t in 1:n) {
#       A[[t]] = A_mat
#     }
#   }
#
#   x_tt_1 = NA
#   P_tt_1 = NA
#   y_tt_1 = NA
#   S_t = NA
#   objective = 0
#
#   if (!return_obj) {
#     filtered_states = matrix(0, nrow = dim_state, ncol = n)
#     filtered_observations = matrix(0, nrow = dim_obs, ncol = n)
#     predicted_states = matrix(0, nrow = dim_state, ncol = n)
#     predicted_observations = matrix(0, nrow = dim_obs, ncol = n)
#     predicted_observations_var = list()
#     mahalanobis_residuals = NA
#   }
#
#   for (t in 1:n) {
#     x_tt_1 = Phi %*% x_tt
#     P_tt_1 = Phi %*% P_tt %*% t(Phi) + Sigma_w
#     y_tt_1 = A[[t]] %*% x_tt_1
#     S_t = A[[t]] %*% P_tt_1 %*% t(A[[t]]) + Sigma_v
#     inv_S_t = solve(S_t)
#     if (sum(abs(gamma[,t])) == 0) {
#       K_t = P_tt_1 %*% t(A[[t]]) %*% inv_S_t
#       x_tt = x_tt_1 + K_t %*% (y[,t] - y_tt_1)
#       P_tt = P_tt_1 - K_t %*% A[[t]] %*% P_tt_1
#     } else {
#       x_tt = x_tt_1
#       P_tt = P_tt_1
#     }
#     if (return_obj) {
#       objective = objective + 1/(2*n) * ((sum(abs(gamma[,t])) == 0) * log(det(S_t)) + t(y[,t] - y_tt_1 - gamma[,t]) %*% inv_S_t %*% (y[,t] - y_tt_1 - gamma[,t]))
#     } else {
#       filtered_states[,t] = x_tt
#       filtered_observations[,t] = A[[t]] %*% x_tt
#       predicted_states[,t] = x_tt_1
#       predicted_observations[,t] = y_tt_1
#       predicted_observations_var[[t]] = S_t
#       mahalanobis_residuals[t] = drop(sqrt(t(y[,t] - y_tt_1) %*% inv_S_t %*% (y[,t] - y_tt_1)))
#     }
#   }
#
#   if (return_obj) {
#     return(objective)
#   } else {
#     return(list(
#       "filtered_states" = filtered_states,
#       "filtered_observations" = filtered_observations,
#       "predicted_states" = predicted_states,
#       "predicted_observations" = predicted_observations,
#       "predicted_observations_var" = predicted_observations_var,
#       "mahalanobis_residuals" = mahalanobis_residuals
#     ))
#   }
# }





# Original used in simulations
# fn_filter = function(
#     par,
#     y,
#     gamma,
#     build,
#     return_obj = FALSE
#     ) {
#
#   SSM_specs = build(par)
#
#   Phi = SSM_specs$state_transition_matrix
#   Sigma_w = SSM_specs$state_noise_var
#   A = SSM_specs$observation_matrix
#   Sigma_v = SSM_specs$observation_noise_var
#   x_tt = SSM_specs$init_state_mean
#   P_tt = SSM_specs$init_state_var
#
#   n = ncol(y)
#   dim_obs = nrow(y)
#   dim_state = nrow(Phi)
#
#   x_tt_1 = NA
#   P_tt_1 = NA
#   y_tt_1 = NA
#   S_t = NA
#   objective = 0
#
#   if (!return_obj) {
#     filtered_states = matrix(0, nrow = dim_state, ncol = n)
#     filtered_observations = matrix(0, nrow = dim_obs, ncol = n)
#     predicted_states = matrix(0, nrow = dim_state, ncol = n)
#     predicted_observations = matrix(0, nrow = dim_obs, ncol = n)
#     predicted_observations_var = list()
#     mahalanobis_residuals = NA
#   }
#
#   for (t in 1:n) {
#     x_tt_1 = Phi %*% x_tt
#     P_tt_1 = Phi %*% P_tt %*% t(Phi) + Sigma_w
#     y_tt_1 = A %*% x_tt_1
#     S_t = A %*% P_tt_1 %*% t(A) + Sigma_v
#     inv_S_t = solve(S_t)
#     if (sum(abs(gamma[,t])) == 0) {
#       K_t = P_tt_1 %*% t(A) %*% inv_S_t
#       x_tt = x_tt_1 + K_t %*% (y[,t] - y_tt_1)
#       P_tt = P_tt_1 - K_t %*% A %*% P_tt_1
#     } else {
#       x_tt = x_tt_1
#       P_tt = P_tt_1
#     }
#     if (return_obj) {
#       objective = objective + 1/(2*n) * ((sum(abs(gamma[,t])) == 0) * log(det(S_t)) + t(y[,t] - y_tt_1 - gamma[,t]) %*% inv_S_t %*% (y[,t] - y_tt_1 - gamma[,t]))
#     } else {
#       filtered_states[,t] = x_tt
#       filtered_observations[,t] = A %*% x_tt
#       predicted_states[,t] = x_tt_1
#       predicted_observations[,t] = y_tt_1
#       predicted_observations_var[[t]] = S_t
#       mahalanobis_residuals[t] = drop(sqrt(t(y[,t] - y_tt_1) %*% inv_S_t %*% (y[,t] - y_tt_1)))
#     }
#   }
#
#   if (return_obj) {
#     return(objective)
#   } else {
#     return(list(
#       "filtered_states" = filtered_states,
#       "filtered_observations" = filtered_observations,
#       "predicted_states" = predicted_states,
#       "predicted_observations" = predicted_observations,
#       "predicted_observations_var" = predicted_observations_var,
#       "mahalanobis_residuals" = mahalanobis_residuals
#     ))
#   }
# }




