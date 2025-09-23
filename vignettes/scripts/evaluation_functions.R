# evaluation functions

get_MSFE = function(model, method, y_oos) {

  x0_oos = c(y_oos[1,], y_oos[1,])
  build_oos = function(par) {

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
  build_oos_huber_trimmed = function(par) {

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
      init_state_mean = x0_oos + c(0, 0, 0.01, 0.01),  # only adjustment
      init_state_var = P0)
  }

  if (method == "classical_model") {
    output = oos_filter(y_oos, model, build_oos)
  } else if (method == "roams_fut_model") {
    output = oos_filter(y_oos, model, build_oos,
                        multiplier = 2)
  } else if (method == "roams_kalman_model") {
    output = oos_filter(y_oos, model, build_oos,
                        multiplier = 1, threshold = Inf)
  } else if (method == "oracle_model") {
    output = oos_filter(y_oos, model, build_oos)
  } else if (method == "huber_model") {
    output = oos_filter(y_oos, model, build_oos_huber_trimmed)
  } else if (method == "trimmed_model") {
    output = oos_filter(y_oos, model, build_oos_huber_trimmed)
  }

  pred = output$predicted_observations
  SFEs = rowSums((pred - y_oos)^2)
  MSFE = mean(SFEs)
  return(MSFE)
}
get_MSFE_contaminated = function(model,
                                 method,
                                 y_oos_contaminated,
                                 y_oos,
                                 true_oos_outlier_locs) {
  x0_oos = c(y_oos_contaminated[1,], y_oos_contaminated[1,])
  build_oos = function(par) {

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
  build_oos_huber_trimmed = function(par) {

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
      init_state_mean = x0_oos + c(0, 0, 0.01, 0.01),  # only adjustment
      init_state_var = P0)
  }

  if (method == "classical_model") {
    output = oos_filter(y_oos_contaminated, model, build_oos)
  } else if (method == "roams_fut_model") {
    output = oos_filter(y_oos_contaminated, model, build_oos,
                        multiplier = 2)
  } else if (method == "roams_kalman_model") {
    output = oos_filter(y_oos_contaminated, model, build_oos,
                        multiplier = 1, threshold = Inf)
  } else if (method == "oracle_model") {
    y_oos_oracle = y_oos_contaminated
    y_oos_oracle[true_oos_outlier_locs,] = NA
    output = oos_filter(y_oos_oracle, model, build_oos)
  } else if (method == "huber_model") {
    output = oos_filter(y_oos_contaminated, model, build_oos_huber_trimmed)
  } else if (method == "trimmed_model") {
    output = oos_filter(y_oos_contaminated, model, build_oos_huber_trimmed)
  }

  pred = output$predicted_observations
  SFEs = rowSums((pred - y_oos)^2)
  MSFE = mean(SFEs)
  return(MSFE)
}

get_TPR = function(model, outlier_locs) {
  if (!inherits(model, "roams_SSM")) {
    return(NA)
  }
  if (sum(outlier_locs == 1) == 0) {
    return(NA)
  }

  nz = rowSums(model$gamma)
  flagged = ifelse(nz != 0, 1, 0)

  TPR = sum((outlier_locs == 1) & (flagged == 1)) / sum(outlier_locs == 1)
  return(TPR)
}
get_TNR = function(model, outlier_locs) {
  if (!inherits(model, "roams_SSM")) {
    return(NA)
  }

  nz = rowSums(model$gamma)
  flagged = ifelse(nz != 0, 1, 0)

  TNR = sum((outlier_locs == 0) & (flagged == 0)) / sum(outlier_locs == 0)
  return(TNR)
}
get_TPR_multilevel = function(model, outlier_levels, level) {
  if (class(model) != "roams_SSM") {
    return(NA)
  }
  if (is.null(outlier_levels)) {
    return(NA)
  }

  nz = rowSums(model$gamma)
  flagged = ifelse(nz != 0, 1, 0)

  TPR = sum((outlier_levels == level) & (flagged == 1)) / sum(outlier_levels == level)
  return(TPR)
}

get_par_dist_phi = function(model, true_par) {
  par_dist = sqrt(sum((model$par[1] - true_par[1])^2))
  return(par_dist)
}
get_par_dist_obs_var = function(model, true_par) {
  par_dist = sqrt(sum((model$par[4:5] - true_par[4:5])^2))
  return(par_dist)
}
get_par_dist_state_var = function(model, true_par) {
  par_dist = sqrt(sum((model$par[2:3] - true_par[2:3])^2))
  return(par_dist)
}
