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
oracle_SSM = function(
  y,
  init_par,
  build,
  outlier_locs,
  B = 20,
  lower = NA,
  upper = NA
  ) {

  if (is.na(lower)[1]) {lower = rep(-Inf, length(init_par))}
  if (is.na(upper)[1]) {upper = rep(Inf, length(init_par))}

  n = ncol(y)
  dim_obs = nrow(y)
  par = init_par
  gamma = matrix(0, nrow = dim_obs, ncol = n)
  r = NA
  for (j in 1:B) {
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

    par = res$par
    filter_output = fn_filter(res$par, y, gamma, build)
    r = y - filter_output$predicted_observations
    gamma_old = gamma
    gamma = matrix(0, nrow = dim_obs, ncol = n)
    gamma[,which(outlier_locs == 1)] = r[,which(outlier_locs == 1)]
    gap = max(abs(gamma - gamma_old))

    nz = sum(colSums(abs(gamma_old)) != 0)
    prop_outlying = nz / n

    if (gap < 1e-4) {
      break
    }
    # new termination criterion
    if (prop_outlying >= 0.5) {
      break
    }
  }

  IPOD_output = list(
    "par" = res$par,
    "gamma" = gamma_old,
    "iterations" = j
  )

  model = c(IPOD_output, filter_output)
  class(model) = "oracle_SSM"
  return(model)
}

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
no_gamma_oracle_SSM = function(
    y,
    init_par,
    build,
    outlier_locs,
    lower = NA,
    upper = NA
) {

  if (is.na(lower)[1]) {lower = rep(-Inf, length(init_par))}
  if (is.na(upper)[1]) {upper = rep(Inf, length(init_par))}

  res = stats::optim(
    par = init_par,
    fn = no_gamma_oracle_filter,
    y = y,
    outlier_locs = outlier_locs,
    build = build,
    return_obj = TRUE,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper
  )

  filter_output = no_gamma_oracle_filter(res$par, y, outlier_locs, build)
  optim_output = list(
    "par" = res$par
  )

  model = c(optim_output, filter_output)
  class(model) = "no_gamma_oracle_SSM"
  return(model)
}

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
delete_outliers_oracle_SSM = function(
    y,
    init_par,
    build,
    outlier_locs,
    lower = NA,
    upper = NA
) {

  if (is.na(lower)[1]) {lower = rep(-Inf, length(init_par))}
  if (is.na(upper)[1]) {upper = rep(Inf, length(init_par))}

  if (sum(outlier_locs) > 0) {
    y_deleted_outliers = y[,-which(outlier_locs == 1)]
  } else {
    y_deleted_outliers = y
  }

  n = ncol(y_deleted_outliers)
  dim_obs = nrow(y_deleted_outliers)
  gamma = matrix(0, nrow = dim_obs, ncol = n)

  res = stats::optim(
    par = init_par,
    fn = fn_filter,
    y = y_deleted_outliers,
    gamma = gamma,
    build = build,
    return_obj = TRUE,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper
  )

  filter_output = no_gamma_oracle_filter(res$par, y, outlier_locs, build)
  optim_output = list(
    "par" = res$par
  )

  model = c(optim_output, filter_output)
  class(model) = "delete_outliers_oracle_SSM"
  return(model)
}

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
classical_SSM = function(
    y,
    init_par,
    build,
    lower = NA,
    upper = NA
) {

  if (is.na(lower)[1]) {lower = rep(-Inf, length(init_par))}
  if (is.na(upper)[1]) {upper = rep(Inf, length(init_par))}

  n = ncol(y)
  dim_obs = nrow(y)
  gamma = matrix(0, nrow = dim_obs, ncol = n)
  res = stats::optim(
    par = init_par,
    fn = fn_filter,
    y = y,
    gamma = gamma,
    build = build,
    return_obj = TRUE,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper
  )

  optim_output = list("par" = res$par)

  filter_output = fn_filter(res$par,
                            y,
                            gamma,
                            build)

  model = c(optim_output, filter_output)
  class(model) = "classical_SSM"
  return(model)
}

#' #' A Cat Function
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
huber_robust_SSM = function(
    y,
    init_par,
    build,
    lower = NA,
    upper = NA
    ) {

  if (is.na(lower)[1]) {lower = rep(-Inf, length(init_par))}
  if (is.na(upper)[1]) {upper = rep(Inf, length(init_par))}

  res = stats::optim(
    par = init_par,
    fn = ruben_filter,
    y = y,
    build = build,
    return_obj = TRUE,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper
  )

  optim_output = list("par" = res$par)

  filter_output = ruben_filter(res$par,
                               y,
                               build)

  model = c(optim_output, filter_output)
  class(model) = "huber_robust_SSM"
  return(model)
}

psi_huber = function(x, k = 2) {ifelse(abs(x) < k, x, sign(x)*k)}  # Ruben uses k = 2
rho_huber_mv = function(x, k = sqrt(qchisq(0.95, df = length(x)))) {
  # df parameter in qchi is the dimension of the observation vector
  norm_x = norm(as.matrix(x), type = "2")
  if (norm_x < k) {
    1/2 * norm_x^2
  } else {
    k*norm_x - 1/2 * k^2
  }
}

ruben_filter = function(
    par,
    y,
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

  d = dim_obs  # dimension of the observation vector
  k = sqrt(qchisq(0.95, d))
  c_H = d / (d*pchisq(k^2, df = d+2) + 2*k*sqrt(2)*gamma((d+1)/2)/gamma(d/2)*(1-pchisq(k^2, df = d+1)) - k^2*(1-pchisq(k^2, df = d)))

  x_tt_1 = NA
  P_tt_1 = NA
  y_tt_1 = NA
  S_t = NA
  objective = 0

  if (!return_obj) {
    filtered_states = matrix(0, nrow = dim_state, ncol = n)
    filtered_observations = matrix(0, nrow = dim_obs, ncol = n)
    predicted_states = matrix(0, nrow = dim_state, ncol = n)
    predicted_observations = matrix(0, nrow = dim_obs, ncol = n)
    predicted_observations_var = list()
    mahalanobis_residuals = NA
  }

  for (t in 1:n) {
    x_tt_1 = Phi %*% x_tt
    P_tt_1 = Phi %*% P_tt %*% t(Phi) + Sigma_w
    y_tt_1 = A %*% x_tt_1

    sqrt_Sigma_v = expm::sqrtm(Sigma_v)
    inv_sqrt_Sigma_v = solve(sqrt_Sigma_v)
    W_elements = drop(psi_huber(inv_sqrt_Sigma_v %*% (y[,t] - y_tt_1))) / drop(inv_sqrt_Sigma_v %*% (y[,t] - y_tt_1))
    inv_W = diag(1/W_elements)
    S_t = A %*% P_tt_1 %*% t(A) + sqrt_Sigma_v %*% inv_W %*% sqrt_Sigma_v
    inv_S_t = solve(S_t)

    if (is.na(y[1,t])) {
      x_tt = x_tt_1
      P_tt = P_tt_1
    } else {
      K_t = P_tt_1 %*% t(A) %*% inv_S_t
      x_tt = x_tt_1 + K_t %*% (y[,t] - y_tt_1)
      P_tt = P_tt_1 - K_t %*% A %*% P_tt_1
    }

    if (return_obj) {
      objective = objective + 1/(2*n) * log(det(S_t)) + c_H/n * rho_huber_mv(expm::sqrtm(inv_S_t) %*% (y[,t] - y_tt_1))
    } else {
      filtered_states[,t] = x_tt
      filtered_observations[,t] = A %*% x_tt
      predicted_states[,t] = x_tt_1
      predicted_observations[,t] = y_tt_1
      predicted_observations_var[[t]] = S_t
      mahalanobis_residuals[t] = drop(sqrt(t(y[,t] - y_tt_1) %*% inv_S_t %*% (y[,t] - y_tt_1)))
    }
  }

  if (return_obj) {
    return(objective)
  } else {
    return(list(
      "filtered_states" = filtered_states,
      "filtered_observations" = filtered_observations,
      "predicted_states" = predicted_states,
      "predicted_observations" = predicted_observations,
      "predicted_observations_var" = predicted_observations_var,
      "mahalanobis_residuals" = mahalanobis_residuals
    ))
  }
}

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
gamma_start_SSM = function(
    y,
    lambda,
    init_par,
    build,
    gamma = array(0, dim(y)),
    B = 20,
    lower = NA,
    upper = NA
) {

  if (is.na(lower)[1]) {lower = rep(-Inf, length(init_par))}
  if (is.na(upper)[1]) {upper = rep(Inf, length(init_par))}

  n = ncol(y)
  dim_obs = nrow(y)
  par = init_par
  r = NA
  for (j in 1:B) {
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

    # if ((sum(res$par == lower) + sum(res$par == upper)) == 0) {
    #   par = res$par
    # } else {
    #   par = init_par
    # }

    filter_output = fn_filter(res$par, y, gamma, build)
    r = y - filter_output$predicted_observations
    gamma_old = gamma
    gamma = matrix(0, nrow = dim_obs, ncol = n)
    gamma[,filter_output$mahalanobis_residuals > lambda] = r[,filter_output$mahalanobis_residuals > lambda]
    gap = max(abs(gamma - gamma_old))

    nz = sum(colSums(abs(gamma_old)) != 0)
    prop_outlying = nz / n

    if (gap < 1e-4) {
      break
    }
    # new termination criterion
    if (prop_outlying >= 0.5) {
      break
    }
  }

  p = length(init_par)
  RSS = sum((r - gamma_old)^2)
  BIC = (n-p)*log(RSS/(n-p)) + (nz+1)*(log(n-p) + 1)
  #negloglik = n*fn_filter(model$par, gamma = gamma_old, y = y, return_obj = TRUE)

  return(list(
    "lambda" = lambda,
    "par" = res$par,
    "prop_outlying" = prop_outlying,
    "BIC" = BIC,
    "RSS" = RSS,
    "gamma" = gamma_old,
    "iterations" = j,
    "r" = r,
    "predicted_observations_var" = filter_output$predicted_observations_var
  ))
}

no_gamma_oracle_filter = function(
    par,
    y,
    outlier_locs,
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

  if (!return_obj) {
    filtered_states = matrix(0, nrow = dim_state, ncol = n)
    filtered_observations = matrix(0, nrow = dim_obs, ncol = n)
    predicted_states = matrix(0, nrow = dim_state, ncol = n)
    predicted_observations = matrix(0, nrow = dim_obs, ncol = n)
    predicted_observations_var = list()
    mahalanobis_residuals = NA
  }

  for (t in 1:n) {
    x_tt_1 = Phi %*% x_tt
    P_tt_1 = Phi %*% P_tt %*% t(Phi) + Sigma_w
    y_tt_1 = A %*% x_tt_1
    S_t = A %*% P_tt_1 %*% t(A) + Sigma_v
    inv_S_t = solve(S_t)
    if (outlier_locs[t] == 0) {
      K_t = P_tt_1 %*% t(A) %*% inv_S_t
      x_tt = x_tt_1 + K_t %*% (y[,t] - y_tt_1)
      P_tt = P_tt_1 - K_t %*% A %*% P_tt_1
    } else {
      x_tt = x_tt_1
      P_tt = P_tt_1
    }
    if (return_obj) {
      objective = objective + 1/(2*n) * (outlier_locs[t] == 0) * (log(det(S_t)) + t(y[,t] - y_tt_1) %*% inv_S_t %*% (y[,t] - y_tt_1))
    } else {
      filtered_states[,t] = x_tt
      filtered_observations[,t] = A %*% x_tt
      predicted_states[,t] = x_tt_1
      predicted_observations[,t] = y_tt_1
      predicted_observations_var[[t]] = S_t
      mahalanobis_residuals[t] = drop(sqrt(t(y[,t] - y_tt_1) %*% inv_S_t %*% (y[,t] - y_tt_1)))
    }
  }

  if (return_obj) {
    return(objective)
  } else {
    return(list(
      "filtered_states" = filtered_states,
      "filtered_observations" = filtered_observations,
      "predicted_states" = predicted_states,
      "predicted_observations" = predicted_observations,
      "predicted_observations_var" = predicted_observations_var,
      "mahalanobis_residuals" = mahalanobis_residuals
    ))
  }
}

