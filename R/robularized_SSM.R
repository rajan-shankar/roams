#' Robust Regularized Fitting of State Space Models
#'
#' Fits a robust state space model to multivariate time series data using iterative parameter estimation and outlier detection. This procedure is inspired by the Iterative Procedure for Outlier Detection (IPOD) algorithm of She and Owen (2011) and is applied over a sequence of regularization parameters (\eqn{\lambda}'s), identifying outliers via Mahalanobis residuals and re-fitting the model iteratively.
#'
#' @param y A numeric matrix of observations, with each row corresponding to a time point.
#' @param init_par A numeric vector of initial parameter values for optimization.
#' @param build A function that accepts a parameter vector and returns a \code{dlm} model (as used in \code{dlm::dlmMLE()}). The \code{specify_SSM} function can be used to create this \code{build} function.
#' @param num_lambdas Integer. The number of \eqn{\lambda} values to evaluate. Ignored if \code{custom_lambdas} is specified. Default is 20.
#' @param custom_lambdas Optional numeric vector. If supplied, these are the exact \eqn{\lambda} values used for model fitting. If not provided or set to \code{NA}, then \code{num_lambdas} \eqn{\lambda}'s are automatically chosen.
#' @param cores Integer. Number of CPU cores to use for parallel processing. Default is 1 (sequential execution).
#' @param B Integer. Maximum number of IPOD iterations per \eqn{\lambda}. Default is 50.
#' @param lower Optional numeric vector of lower bounds for optimization. If \code{NA}, defaults to \code{-Inf} for all parameters. Must be of same length as \code{init_par}.
#' @param upper Optional numeric vector of upper bounds for optimization. If \code{NA}, defaults to \code{Inf} for all parameters. Must be of same length as \code{init_par}.
#' @param control A named list of control options to pass to \code{optim} via \code{dlm::dlmMLE()}. Default is \code{list(parscale = init_par)}, which can help the optimizer if parameters are on vastly different scales.
#'
#' @return If more than one \eqn{\lambda} values are used, returns an object of class \code{robularized_SSM_list} — a list containing a \code{robularized_SSM} model for each \eqn{\lambda}. If only one \eqn{\lambda} value is used (i.e. \code{custom_lambdas} is manually specified as a single value), returns a single \code{robularized_SSM} object.
#'
#' Each \code{robularized_SSM} object includes:
#' \itemize{
#'   \item \code{lambda} - The \eqn{\lambda} value used.
#'   \item \code{prop_outlying} - Proportion of non-missing time points identified as outliers.
#'   \item \code{BIC} - Bayesian Information Criterion of the final model.
#'   \item \code{loglik} - Log-likelihood of the fitted model.
#'   \item \code{RSS} - Residual sum of squares.
#'   \item \code{gamma} - Matrix of estimated outlier adjustments.
#'   \item \code{iterations} - Number of IPOD iterations performed.
#'   \item Optimization output from \code{dlm::dlmMLE()} from the final IPOD iteration.
#'   \item \code{y} - The original data matrix.
#'   \item \code{build} - The original build function used to specify the model.
#' }
#'
#' @details
#' The IPOD procedure alternates between estimating model parameters via maximum likelihood and identifying outlying observations based on Mahalanobis residuals. For each iteration:
#' \enumerate{
#'   \item A \code{dlm} model is fit using \code{dlm::dlmMLE()}.
#'   \item Mahalanobis residuals are computed.
#'   \item Observations with residuals above the current \eqn{\lambda} threshold are treated as missing in the next iteration.
#' }
#'
#' The algorithm stops when the change in parameters and outlier estimates is sufficiently small or if too many outliers are detected (more than 50\% of complete observations).
#'
#' @seealso \code{\link[dlm]{dlmMLE}}, \code{\link{best_BIC_model}}, \code{\link{outlier_target_model}}, \code{\link{get_attribute}}, \code{\link{autoplot.robularized_SSM_list}}, \code{\link{attach_insample_info}}, \code{\link{oos_filter}}, \code{\link{specify_SSM}}
#'
#' @references She, Y., & Owen, A. B. (2011). Outlier Detection Using Nonconvex Penalized Regression. *Journal of the American Statistical Association, 106*(494), 626–639. https://doi.org/10.1198/jasa.2011.tm10390
#'
#' @export
robularized_SSM = function(
    y,
    init_par,
    build,
    num_lambdas = 20,
    custom_lambdas = NA,
    cores = 1,
    B = 50,
    lower = NA,
    upper = NA,
    control = list(parscale = init_par)
    ) {

  if (ncol(y) > nrow(y)) {
    warning("Input data has more columns than rows. Did you forget to transpose your data? This function expects each row to represent a time point (i.e., observations in rows).")
  }

  if (any(is.na(custom_lambdas))) {
  # Completely automatic fit

    # Classical fit by using large lambda
    classical = run_IPOD(y = y,
                         lambda = 100,
                         init_par = init_par,
                         build = build,
                         B = B,
                         lower = lower,
                         upper = upper,
                         control = control)

    # Highest lambda is the supremum norm of mahalanobis residuals of classical fit
    highest_lambda = max(dlmInfo(y, y, classical, build)$mahalanobis_residuals)
    lowest_lambda = 2

    # Lambda grid
    lambdas = seq(lowest_lambda,
                  highest_lambda,
                  length.out = num_lambdas)

  } else {
  # Lambdas have been manually specified

    if (length(custom_lambdas) != 1) {
      lambdas = custom_lambdas[order(custom_lambdas)]  # ensure lambdas are in ascending order
    } else {
      lambda = custom_lambdas
      model = run_IPOD(y = y,
                       lambda = lambda,
                       init_par = init_par,
                       build = build,
                       B = B,
                       lower = lower,
                       upper = upper,
                       control = control)

      return(model)
    }
  }


  # Fit models across the grid
  model_list = lambda_grid(y = y,
                           lambdas = lambdas,
                           init_par = init_par,
                           build = build,
                           cores = cores,
                           B = B,
                           lower = lower,
                           upper = upper,
                           control = control)

  return(model_list)
}

lambda_grid = function(
    y,
    lambdas,
    init_par,
    build,
    cores,
    B,
    lower,
    upper,
    control
    ) {

  if (cores == 1) {
    model_list = list()
    for (i in 1:length(lambdas)) {
      model_list[[i]] = run_IPOD(y,
                                 lambdas[i],
                                 init_par,
                                 build,
                                 B,
                                 lower,
                                 upper,
                                 control)
    }
  } else {

    future::plan(future::multisession, workers = cores)

    model_list = furrr::future_map(lambdas, ~ run_IPOD(y,
                                                       .x,
                                                       init_par,
                                                       build,
                                                       B,
                                                       lower,
                                                       upper,
                                                       control))

    future::plan(future::sequential)
  }

  #   cl = parallel::makeCluster(cores)
  #   doParallel::registerDoParallel(cl)
  #                                  #export = list(run_IPOD = run_IPOD))
  #   model_list = foreach::foreach(
  #     i = 1:length(lambdas)) %dopar% {
  #                   model = run_IPOD(y,
  #                                    lambdas[i],
  #                                    init_par,
  #                                    build,
  #                                    B,
  #                                    lower,
  #                                    upper,
  #                                    control)
  #                   return(model)
  #                   }
  #   parallel::stopCluster(cl)
  # }

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
    upper,
    control
) {

  if (is.na(lower)[1]) {lower = rep(-Inf, length(init_par))}
  if (is.na(upper)[1]) {upper = rep(Inf, length(init_par))}

  n = nrow(y)
  n_complete = sum(!is.na(rowSums(y)))
  dim_obs = ncol(y)
  par = init_par
  gamma = matrix(0, nrow = n, ncol = dim_obs)
  adj_y = y
  r = NA
  theta_old = par

  for (j in 1:B) {
    if (j != 1) {theta_old = fit$par}

    fit = dlm::dlmMLE(
      adj_y,
      parm = par,
      build = build,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = control
      )

    # par = init_par


    # If the fit is at edge of parameter space, use the initial parameters for next iteration
    if ((sum(fit$par == lower) + sum(fit$par == upper)) == 0) {
      par = fit$par
    } else {
      par = init_par
    }

    # Update gammas
    info_output = dlmInfo(y, adj_y, fit, build)
    r = y - info_output$predicted_observations
    gamma_old = gamma
    gamma = matrix(0, nrow = n, ncol = dim_obs)
    gamma[info_output$mahalanobis_residuals > lambda,] = r[info_output$mahalanobis_residuals > lambda,]
    adj_y = y
    adj_y[which(rowSums(abs(gamma)) != 0),] = NA
    gap_gamma = max(abs(gamma - gamma_old))
    gap_theta = max(abs(fit$par - theta_old))

    nz = sum(rowSums(abs(gamma_old)) != 0)
    prop_outlying = nz / n_complete

    if ((gap_gamma < 1e-4) && (gap_theta < 1e-4)) {
      break
    }

    # Too many outliers detected will cause instability
    if (prop_outlying >= 0.5) {
      break
    }
  }

  p = length(init_par)
  RSS = sum((r - gamma_old)^2, na.rm = TRUE)
  loglik = -dlm::dlmLL(adj_y, mod = build(fit$par))
  BIC = nz*log(n_complete) - 2*loglik

  model = c(
    list(
      "lambda" = lambda,
      "prop_outlying" = prop_outlying,
      "BIC" = BIC,
      "loglik" = loglik,
      "RSS" = RSS,
      "gamma" = gamma_old,
      "iterations" = j
      ),
    fit,  # output from dlmMLE (which is just output from optim)
    list(y = y),
    list(build = build)
    )

  class(model) = "robularized_SSM"
  return(model)
}

IPOD_oos_robust_filter = function(y, par, build, threshold, multiplier = 1) {

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

  x_tt_1 = NA
  P_tt_1 = NA
  y_tt_1 = NA
  S_t = NA

  filtered_states = matrix(0, nrow = n, ncol = dim_state)
  filtered_observations = matrix(0, nrow = n, ncol = dim_obs)
  predicted_states = matrix(0, nrow = n, ncol = dim_state)
  predicted_observations = matrix(0, nrow = n, ncol = dim_obs)
  filtered_states_var = list()
  predicted_states_var = list()
  predicted_observations_var = list()
  mahalanobis_residuals = NA
  outliers_flagged = rep(0, n)

  for (t in 1:n) {
    x_tt_1 = Phi %*% x_tt
    P_tt_1 = Phi %*% P_tt %*% t(Phi) + Sigma_w
    y_tt_1 = A %*% x_tt_1
    S_t = A %*% P_tt_1 %*% t(A) + Sigma_v
    inv_S_t = solve(S_t)

    if (any(is.na(y[t,]))) {
      mahalanobis_residuals[t] = 0
      x_tt = x_tt_1
      P_tt = multiplier*P_tt_1
    } else {

      mahalanobis_residuals[t] = drop(sqrt(t(y[t,] - y_tt_1) %*% inv_S_t %*% (y[t,] - y_tt_1)))

      if (mahalanobis_residuals[t] <= threshold) {
        K_t = P_tt_1 %*% t(A) %*% inv_S_t
        x_tt = x_tt_1 + K_t %*% (y[t,] - y_tt_1)
        P_tt = P_tt_1 - K_t %*% A %*% P_tt_1
      } else {
        x_tt = x_tt_1
        P_tt = multiplier*P_tt_1
        outliers_flagged[t] = 1
      }
    }

    filtered_states[t,] = x_tt
    filtered_observations[t,] = A %*% x_tt
    predicted_states[t,] = x_tt_1
    predicted_observations[t,] = y_tt_1
    filtered_states_var[[t]] = P_tt
    predicted_states_var[[t]] = P_tt_1
    predicted_observations_var[[t]] = S_t
  }

  return(list(
    "filtered_states" = filtered_states,
    "filtered_observations" = filtered_observations,
    "predicted_states" = predicted_states,
    "predicted_observations" = predicted_observations,
    "filtered_states_var" = filtered_states_var,
    "predicted_states_var" = predicted_states_var,
    "predicted_observations_var" = predicted_observations_var,
    "mahalanobis_residuals" = mahalanobis_residuals,
    "outliers_flagged" = outliers_flagged
  ))
}

# A simpler version of attach_insample_info for internal use only.
# Mainly used for getting residuals and predictions from robularized models.
dlmInfo = function(y, adj_y, model, build) {

  filter_output = dlm::dlmFilter(adj_y, mod = build(model$par))
  smoother_output = dlm::dlmSmooth(filter_output)
  A = filter_output$mod$FF

  S = purrr::map(dlm::dlmSvd2var(filter_output$U.R, filter_output$D.R),
          ~ A %*% . %*% t(A) + filter_output$mod$V)
  inv_S = purrr::map(S, ~ solve(.))
  mahalanobis_residuals = purrr::map2_dbl(
    apply(y - filter_output$f, 1, c, simplify = FALSE),
    inv_S,
    ~ drop(t(.x) %*% .y %*% .x)) %>% sqrt()

  mahalanobis_residuals = ifelse(is.na(mahalanobis_residuals), 0, mahalanobis_residuals)

  return(list(
    smoothed_observations = (smoother_output$s %*% t(A))[2:(nrow(y) + 1),],
    filtered_observations = (filter_output$m %*% t(A))[2:(nrow(y) + 1),],
    predicted_observations = filter_output$f,
    mahalanobis_residuals = mahalanobis_residuals
  ))

}

#' Attach In-Sample Information to Fitted State Space Model
#'
#' Attaches detailed in-sample information—such as predicted, filtered, and smoothed states and observations—to a model object fitted using any of the package’s supported SSM estimation methods.
#' These quantities are not stored by default in model objects due to their potentially large memory footprint.
#'
#' @param model A fitted model object of class \code{robularized_SSM}, \code{classical_SSM}, \code{oracle_SSM}, \code{huber_robust_SSM}, or \code{trimmed_robust_SSM}.
#'
#' @return A modified version of the input model object, with an additional class \code{insample_info}, and the following in-sample elements appended:
#' \describe{
#'   \item{\code{filtered_states}}{Filtered state estimates using data up to each time point.}
#'   \item{\code{predicted_states}}{One-step-ahead state predictions.}
#'   \item{\code{filtered_observations}}{Expected observations given data up to each time point.}
#'   \item{\code{predicted_observations}}{One-step-ahead forecasts of observations.}
#'   \item{\code{filtered_states_var}}{List of filtered state variance matrices.}
#'   \item{\code{predicted_states_var}}{List of one-step-ahead state prediction variances.}
#'   \item{\code{predicted_observations_var}}{List of one-step-ahead observation forecast variances.}
#'   \item{\code{mahalanobis_residuals}}{Vector of Mahalanobis distances of residuals from predicted observations.}
#' }
#'
#' For models of class \code{robularized_SSM}, \code{classical_SSM}, or \code{oracle_SSM}, the following additional elements are also attached:
#' \describe{
#'   \item{\code{smoothed_states}}{Posterior means of hidden states using all data.}
#'   \item{\code{smoothed_observations}}{Posterior mean of the observed series based on smoothed states.}
#'   \item{\code{smoothed_states_var}}{List of smoothed state variance matrices.}
#' }
#'
#' @details
#' The attached outputs enable richer diagnostics, outlier inspection, and plotting.
#' For \code{huber_robust_SSM} and \code{trimmed_robust_SSM} models, in-sample information is computed using a custom robust filtering function, and smoothed quantities (\code{smoothed_states}, \code{smoothed_observations}, and \code{smoothed_states_var}) are \strong{not available}.
#' This function should only be applied once to a model object.
#'
#' @seealso \code{\link{oos_filter}}
#'
#' @import tidyverse
#'
#' @export
attach_insample_info = function(model) {

  if (inherits(model, "insample_info")) {
    stop("Model already has in-sample information attached.")
  }

  y = model$y

  if (inherits(model, "huber_robust_SSM")) {
    insample_info = ruben_filter(model$par, y, model$build, obj_type = "huber")
    output = c(model, insample_info)
    class(output) = c("huber_robust_SSM", "insample_info")
    return(output)
  } else if (inherits(model, "trimmed_robust_SSM")) {
    insample_info = ruben_filter(model$par, y, model$build, obj_type = "trimmed", alpha = model$alpha)
    output = c(model, insample_info)
    class(output) = c("trimmed_robust_SSM", "insample_info")
    return(output)

  } else if (inherits(model, "robularized_SSM")){
    adj_y = y
    adj_y[which(rowSums(abs(model$gamma)) != 0),] = NA
  } else if (inherits(model, "classical_SSM")) {
    adj_y = y
  } else if (inherits(model, "oracle_SSM")) {
    adj_y = y
    adj_y[model$outlier_locs != 0,] = NA
  } else {
    stop("Invalid model class. Expected 'robularized_SSM' or 'classical_SSM' or 'oracle_SSM' or 'huber_robust_SSM' or 'trimmed_robust_SSM'.")
  }

  filter_output = dlm::dlmFilter(adj_y, mod = model$build(model$par))
  smoother_output = dlm::dlmSmooth(filter_output)
  A = filter_output$mod$FF

  P_tt_1 = dlm::dlmSvd2var(filter_output$U.R, filter_output$D.R)
  S = purrr::map(P_tt_1, ~ A %*% . %*% t(A) + filter_output$mod$V)
  P_tt = dlm::dlmSvd2var(filter_output$U.C, filter_output$D.C)
  P_tt_smooth = dlm::dlmSvd2var(smoother_output$U.S, smoother_output$D.S)

  inv_S = purrr::map(S, ~ solve(.))
  mahalanobis_residuals = purrr::map2_dbl(
    apply(y - filter_output$f, 1, c, simplify = FALSE),
    inv_S,
    ~ drop(t(.x) %*% .y %*% .x)) %>% sqrt()

  mahalanobis_residuals = ifelse(is.na(mahalanobis_residuals), 0, mahalanobis_residuals)

  output = c(model,
    list(
    smoothed_states  = smoother_output$s[2:(nrow(y) + 1),],
    filtered_states  = filter_output$m[2:(nrow(y) + 1),],
    predicted_states = filter_output$a,
    smoothed_observations  = (smoother_output$s %*% t(A))[2:(nrow(y) + 1),],
    filtered_observations  = (filter_output$m %*% t(A))[2:(nrow(y) + 1),],
    predicted_observations = filter_output$f,
    smoothed_states_var = P_tt_smooth,
    filtered_states_var = P_tt,
    predicted_states_var = P_tt_1,
    predicted_observations_var = S,
    mahalanobis_residuals = mahalanobis_residuals
  ))

  if (inherits(model, "robularized_SSM")) {
    class(output) = c("robularized_SSM", "insample_info")
  } else if (inherits(model, "classical_SSM")) {
    class(output) = c("classical_SSM", "insample_info")
  } else if (inherits(model, "oracle_SSM")) {
    class(output) = c("oracle_SSM", "insample_info")
  }

  return(output)
}

#' Compute Out-of-Sample Inference for Fitted State Space Model
#'
#' Applies the fitted model parameters to a user-supplied out-of-sample dataset to compute predicted and filtered states and observations. Robust and classical inference procedures are supported depending on the class of the input model.
#'
#' @param y_oos A numeric matrix containing out-of-sample observations. Each row corresponds to a time point.
#' @param model A fitted model object of class \code{robularized_SSM}, \code{classical_SSM}, \code{oracle_SSM}, \code{huber_robust_SSM}, or \code{trimmed_robust_SSM}.
#' @param build A function that maps a numeric parameter vector to a corresponding \code{dlm} model object. The \code{specify_SSM} function can be used to create this \code{build} function.
#' @param outlier_locs A logical or binary vector of the same length as \code{nrow(y)}, indicating time points to be treated as missing (i.e., time points that are known to be outliers). Used only with \code{oracle_SSM} models.
#' @param threshold Mahalanobis distance threshold for identifying outliers in \code{robularized_SSM} models. Default is \code{sqrt(qchisq(0.99, ncol(y)))}.
#' @param multiplier Multiplier for how quickly the filter grows its filtered state variance (uncertainty) after detecting an outlier in \code{robularized_SSM} models. Default is \code{2}.
#'
#' @return A named list containing out-of-sample inference results:
#' \describe{
#'   \item{\code{filtered_states}}{Filtered state estimates using out-of-sample data.}
#'   \item{\code{predicted_states}}{One-step-ahead state predictions.}
#'   \item{\code{filtered_observations}}{Expected observations given past out-of-sample data.}
#'   \item{\code{predicted_observations}}{One-step-ahead forecasts of observations.}
#'   \item{\code{filtered_states_var}}{List of filtered state variance matrices.}
#'   \item{\code{predicted_states_var}}{List of one-step-ahead state prediction variances.}
#'   \item{\code{predicted_observations_var}}{List of one-step-ahead observation forecast variances.}
#'   \item{\code{mahalanobis_residuals}}{Vector of Mahalanobis distances of residuals from predicted observations.}
#' }
#'
#' @details
#' The function reuses the model's fitted parameters to generate inference on new data \code{y_oos}. Robust variants use appropriate robust filters, while the classical and oracle models use standard Kalman filtering. For \code{oracle_SSM} models, observations flagged in \code{outlier_locs} are treated as missing during filtering.
#'
#' @seealso \code{\link{attach_insample_info}}, \code{\link{specify_SSM}}
#'
#' @import tidyverse
#'
#' @export
oos_filter = function(y_oos, model, build,
                      outlier_locs = rep(0, nrow(y_oos)),
                      threshold = sqrt(qchisq(0.99, ncol(y_oos))),
                      multiplier = 2) {

  y = y_oos

  if (inherits(model, "huber_robust_SSM")) {
    oos_output = ruben_filter(model$par, y, build, obj_type = "huber")
    return(oos_output)
  } else if (inherits(model, "trimmed_robust_SSM")) {
    oos_output = ruben_filter(model$par, y, build, obj_type = "trimmed")
    return(oos_output)

  } else if (inherits(model, "robularized_SSM")){

    # Robust threshold filter
    oos_output = IPOD_oos_robust_filter(
      y = y,
      par = model$par,
      build = build,
      threshold = threshold,
      multiplier = multiplier
    )

    return(oos_output)

  } else if (inherits(model, "classical_SSM")) {
    adj_y = y
  } else if (inherits(model, "oracle_SSM")) {
    adj_y = y
    adj_y[outlier_locs != 0,] = NA
  } else {
    stop("Invalid model class. Expected 'robularized_SSM' or 'classical_SSM' or 'oracle_SSM' or 'huber_robust_SSM' or 'trimmed_robust_SSM'.")
  }

  filter_output = dlm::dlmFilter(adj_y, mod = build(model$par))
  A = filter_output$mod$FF

  P_tt_1 = dlm::dlmSvd2var(filter_output$U.R, filter_output$D.R)
  S = purrr::map(P_tt_1, ~ A %*% . %*% t(A) + filter_output$mod$V)
  P_tt = dlm::dlmSvd2var(filter_output$U.C, filter_output$D.C)

  inv_S = purrr::map(S, ~ solve(.))
  mahalanobis_residuals = purrr::map2_dbl(
    apply(y - filter_output$f, 1, c, simplify = FALSE),
    inv_S,
    ~ drop(t(.x) %*% .y %*% .x)) %>% sqrt()

  mahalanobis_residuals = ifelse(is.na(mahalanobis_residuals), 0, mahalanobis_residuals)

  return(list(
             filtered_states  = filter_output$m[2:(nrow(y) + 1),],
             predicted_states = filter_output$a,
             filtered_observations  = (filter_output$m %*% t(A))[2:(nrow(y) + 1),],
             predicted_observations = filter_output$f,
             filtered_states_var = P_tt,
             predicted_states_var = P_tt_1,
             predicted_observations_var = S,
             mahalanobis_residuals = mahalanobis_residuals
           ))

}



