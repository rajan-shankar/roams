#' Oracle State Space Model Fit
#'
#' Fits a state space model by treating a known set of outliers as missing data. This benchmark model assumes prior knowledge of outlier locations and is intended for comparison with automatic outlier detection procedures.
#'
#' @param y A numeric matrix of observations (time points in rows).
#' @param init_par A numeric vector of initial parameter values.
#' @param build A function that returns a \code{dlm} model given a parameter vector. The \code{specify_SSM()} function can be used to create this \code{build} function.
#' @param outlier_locs An integer or logical vector of length equal to the number of time points, indicating locations of known outliers.
#' @param lower Optional numeric vector of lower bounds for parameter estimation. Defaults to \code{-Inf}. Must be of same length as \code{init_par}.
#' @param upper Optional numeric vector of upper bounds for parameter estimation. Defaults to \code{Inf}. Must be of same length as \code{init_par}.
#' @param control Optional list of control parameters passed to \code{optim} via \code{dlm::dlmMLE()}. Default is \code{list(parscale = init_par)}, which can help the optimizer if parameters are on vastly different scales.
#'
#' @return An object of class \code{oracle_SSM} containing the optimization result, the provided outlier locations, the original data, and the original build function.
#'
#' @seealso \code{\link[dlm]{dlmMLE}}, \code{\link{robularized_SSM}}, \code{\link{attach_insample_info}}, \code{\link{oos_filter}}, \code{\link{specify_SSM}}
#'
#' @export
oracle_SSM = function(
    y,
    init_par,
    build,
    outlier_locs,
    lower = NA,
    upper = NA,
    control = list(parscale = init_par)
) {

  if (is.na(lower)[1]) {lower = rep(-Inf, length(init_par))}
  if (is.na(upper)[1]) {upper = rep(Inf, length(init_par))}

  adj_y = y
  adj_y[outlier_locs != 0,] = NA

  optim_output = dlm::dlmMLE(
    adj_y,
    parm = init_par,
    build,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = control
  )

  model = c(optim_output, "outlier_locs" = outlier_locs, list("y" = y), list("build" = build))
  class(model) = "oracle_SSM"
  return(model)
}

#' Classical State Space Model Fit
#'
#' Fits a state space model using classical maximum likelihood estimation with no attempt to detect or account for outliers. This serves as a baseline model for comparison.
#'
#' @param y A numeric matrix of observations (time points in rows).
#' @param init_par A numeric vector of initial parameter values.
#' @param build A function that returns a \code{dlm} model given a parameter vector. The \code{specify_SSM()} function can be used to create this \code{build} function.
#' @param lower Optional numeric vector of lower bounds for parameter estimation. Defaults to \code{-Inf}. Must be of same length as \code{init_par}.
#' @param upper Optional numeric vector of upper bounds for parameter estimation. Defaults to \code{Inf}. Must be of same length as \code{init_par}.
#' @param control Optional list of control parameters passed to \code{optim} via \code{dlm::dlmMLE()}. Default is \code{list(parscale = init_par)}, which can help the optimizer if parameters are on vastly different scales.
#'
#' @return An object of class \code{classical_SSM} containing the optimization result, the original data, and the original build function.
#'
#' @seealso \code{\link[dlm]{dlmMLE}}, \code{\link{robularized_SSM}}, \code{\link{attach_insample_info}}, \code{\link{oos_filter}}, \code{\link{specify_SSM}}
#'
#' @export
classical_SSM = function(
    y,
    init_par,
    build,
    lower = NA,
    upper = NA,
    control = list(parscale = init_par)
) {

  if (is.na(lower)[1]) {lower = rep(-Inf, length(init_par))}
  if (is.na(upper)[1]) {upper = rep(Inf, length(init_par))}

  n = nrow(y)
  dim_obs = ncol(y)

  optim_output = dlm::dlmMLE(
    y,
    parm = init_par,
    build,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = control
  )

  model = c(optim_output, list("y" = y), list("build" = build))
  class(model) = "classical_SSM"
  return(model)
}

#' Huber-Robust State Space Model Fit
#'
#' Fits a robust state space model by minimizing a Huber loss objective as per Crevits and Croux (2018), providing protection against moderate outliers. The predicted observations used in the Huber loss are computed using the Huber robust filter from Cipra and Romera (1997).
#'
#' @param y A numeric matrix of observations (time points in rows).
#' @param init_par A numeric vector of initial parameter values.
#' @param build A function that returns a \code{dlm} model given a parameter vector. The \code{specify_SSM()} function can be used to create this \code{build} function.
#' @param lower Optional numeric vector of lower bounds for parameter estimation. Defaults to \code{-Inf}. Must be of same length as \code{init_par}.
#' @param upper Optional numeric vector of upper bounds for parameter estimation. Defaults to \code{Inf}. Must be of same length as \code{init_par}.
#' @param control Optional list of control parameters passed to \code{optim}. Default is \code{list(parscale = init_par)}, which can help the optimizer if parameters are on vastly different scales.
#'
#' @return An object of class \code{huber_robust_SSM} containing the optimization result, the original data, and the original build function.
#'
#' @seealso \code{\link{trimmed_robust_SSM}}, \code{\link{robularized_SSM}}, \code{\link[stats]{optim}}, \code{\link{attach_insample_info}}, \code{\link{oos_filter}}, \code{\link{specify_SSM}}
#'
#' @references Crevits R. and Croux C. (2018). Robust Estimation of Linear State Space Models. *Communications in Statistics: Simulation and Computation*
#' @references Cipra, T., Romera, R. (1997). Kalman filter with outliers and missing observations. *Test* 6, 379–395. https://doi.org/10.1007/BF02564705
#'
#' @export
huber_robust_SSM = function(
    y,
    init_par,
    build,
    lower = NA,
    upper = NA,
    control = list(parscale = init_par)
    ) {

  if (is.na(lower)[1]) {lower = rep(-Inf, length(init_par))}
  if (is.na(upper)[1]) {upper = rep(Inf, length(init_par))}

  optim_output = stats::optim(
    par = init_par,
    fn = ruben_filter,
    y = y,
    build = build,
    return_obj = TRUE,
    obj_type = "huber",
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = control
  )

  model = c(optim_output, list("y" = y), list("build" = build))
  class(model) = "huber_robust_SSM"
  return(model)
}

#' Trimmed-Robust State Space Model Fit
#'
#' Fits a robust state space model by minimizing a trimmed loss function as per Crevits and Croux (2018). A fixed proportion of the largest residuals are excluded from the objective, providing robustness against extreme outliers. The predicted observations used in the trimmed loss are computed using the Huber robust filter from Cipra and Romera (1997).
#'
#' @param y A numeric matrix of observations (time points in rows).
#' @param init_par A numeric vector of initial parameter values.
#' @param build A function that returns a \code{dlm} model given a parameter vector. The \code{specify_SSM()} function can be used to create this \code{build} function.
#' @param alpha Numeric value in the interval [0, 1) indicating the trimming proportion (i.e., the proportion of data to exclude as outliers).
#' @param lower Optional numeric vector of lower bounds for parameter estimation. Defaults to \code{-Inf}. Must be of same length as \code{init_par}.
#' @param upper Optional numeric vector of upper bounds for parameter estimation. Defaults to \code{Inf}. Must be of same length as \code{init_par}.
#' @param control Optional list of control parameters passed to \code{optim}. Default is \code{list(parscale = init_par)}, which can help the optimizer if parameters are on vastly different scales.
#'
#' @return An object of class \code{trimmed_robust_SSM} containing the optimization result, trimming level \eqn{\alpha}, the original data, and the original build function.
#'
#' @seealso \code{\link{huber_robust_SSM}}, \code{\link{robularized_SSM}}, \code{\link[stats]{optim}}, \code{\link{attach_insample_info}}, \code{\link{oos_filter}}, \code{\link{specify_SSM}}
#'
#' @references Crevits R. and Croux C. (2018). Robust Estimation of Linear State Space Models. *Communications in Statistics: Simulation and Computation*
#' @references Cipra, T., Romera, R. (1997). Kalman filter with outliers and missing observations. *Test* 6, 379–395. https://doi.org/10.1007/BF02564705
#'
#' @export
trimmed_robust_SSM = function(
    y,
    init_par,
    build,
    alpha,
    lower = NA,
    upper = NA,
    control = list(parscale = init_par)
) {

  if (is.na(lower)[1]) {lower = rep(-Inf, length(init_par))}
  if (is.na(upper)[1]) {upper = rep(Inf, length(init_par))}

  optim_output = stats::optim(
    par = init_par,
    fn = ruben_filter,
    y = y,
    build = build,
    return_obj = TRUE,
    obj_type = "trimmed",
    alpha = alpha,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = control
  )

  model = c(optim_output, "alpha" = alpha, list("y" = y), list("build" = build))
  class(model) = "trimmed_robust_SSM"
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
    return_obj = FALSE,
    obj_type = "huber",
    alpha = 0.1
) {

  SSM_specs = build(par)

  Phi = SSM_specs$GG
  Sigma_w = SSM_specs$W
  A = SSM_specs$FF
  Sigma_v = SSM_specs$V
  x_tt = SSM_specs$m0
  P_tt = SSM_specs$C0

  n = nrow(y)
  dim_obs = ncol(y)
  dim_state = nrow(Phi)


  if (return_obj) {
    d = dim_obs  # dimension of the observation vector
    k = sqrt(qchisq(0.95, d))
    if (obj_type == "huber") {
      c_H = d / (d*pchisq(k^2, df = d+2) + 2*k*sqrt(2)*gamma((d+1)/2)/gamma(d/2)*(1-pchisq(k^2, df = d+1)) - k^2*(1-pchisq(k^2, df = d)))
    } else if (obj_type == "trimmed") {
      c_T = 1 / (pchisq(qchisq(1 - alpha, df = d), df = d+2))
      mahalanobis_residuals = NA
      det_S_t = NA
    }
  }

  x_tt_1 = NA
  P_tt_1 = NA
  y_tt_1 = NA
  S_t = NA
  objective = 0

  if (!return_obj) {
    filtered_states = matrix(0, nrow = n, ncol = dim_state)
    filtered_observations = matrix(0, nrow = n, ncol = dim_obs)
    predicted_states = matrix(0, nrow = n, ncol = dim_state)
    predicted_observations = matrix(0, nrow = n, ncol = dim_obs)
    predicted_observations_var = list()
    filtered_states_var = list()
    predicted_states_var = list()
    mahalanobis_residuals = NA
  }

  for (t in 1:n) {
    x_tt_1 = Phi %*% x_tt
    P_tt_1 = Phi %*% P_tt %*% t(Phi) + Sigma_w
    y_tt_1 = A %*% x_tt_1

    sqrt_Sigma_v = expm::sqrtm(Sigma_v)
    inv_sqrt_Sigma_v = solve(sqrt_Sigma_v)
    W_elements = drop(psi_huber(inv_sqrt_Sigma_v %*% (y[t,] - y_tt_1))) / drop(inv_sqrt_Sigma_v %*% (y[t,] - y_tt_1))
    inv_W = diag(1/W_elements)
    S_t = A %*% P_tt_1 %*% t(A) + sqrt_Sigma_v %*% inv_W %*% sqrt_Sigma_v
    inv_S_t = solve(S_t)

    if (any(is.na(y[t,]))) {
      x_tt = x_tt_1
      P_tt = P_tt_1
    } else {
      K_t = P_tt_1 %*% t(A) %*% inv_S_t
      x_tt = x_tt_1 + K_t %*% (y[t,] - y_tt_1)
      P_tt = P_tt_1 - K_t %*% A %*% P_tt_1
    }

    if (return_obj) {
      if (obj_type == "huber") {
        objective = objective + 1/(2*n) * log(det(S_t)) + c_H/n * rho_huber_mv(expm::sqrtm(inv_S_t) %*% (y[t,] - y_tt_1))
      } else if (obj_type == "trimmed") {
        mahalanobis_residuals[t] = drop(sqrt(t(y[t,] - y_tt_1) %*% inv_S_t %*% (y[t,] - y_tt_1)))
        det_S_t[t] = det(S_t)
      }
    } else {
      filtered_states[t,] = x_tt
      filtered_observations[t,] = A %*% x_tt
      predicted_states[t,] = x_tt_1
      predicted_observations[t,] = y_tt_1
      predicted_observations_var[[t]] = S_t
      filtered_states_var[[t]] = P_tt
      predicted_states_var[[t]] = P_tt_1
      mahalanobis_residuals[t] = drop(sqrt(t(y[t,] - y_tt_1) %*% inv_S_t %*% (y[t,] - y_tt_1)))
    }
  }

  # Compute trimmed likelihood
  if (obj_type == "trimmed" & return_obj) {
    keep_set = which(rank(mahalanobis_residuals) < (1 - alpha)*n)
    objective = 1/(2*n*(1-alpha)) * sum(log(det_S_t[keep_set]) + c_T*mahalanobis_residuals[keep_set]^2)
  }

  if (return_obj) {
    return(objective)
  } else {
    return(list(
      filtered_states = filtered_states,
      predicted_states = predicted_states,
      filtered_observations = filtered_observations,
      predicted_observations = predicted_observations,
      filtered_states_var = filtered_states_var,
      predicted_states_var = predicted_states_var,
      predicted_observations_var = predicted_observations_var,
      mahalanobis_residuals = mahalanobis_residuals
    ))
  }
}


