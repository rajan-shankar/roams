#' Robust Outlier-Adjusted Mean-Shift Estimation of State Space Models
#'
#' Fits a robust state space model to multivariate time series data using iterative parameter estimation and outlier detection. This procedure is inspired by the Iterative Procedure for Outlier Detection (IPOD) algorithm of She and Owen (2011) and is applied over a sequence of regularization parameters (\eqn{\lambda}'s), identifying outliers via Mahalanobis residuals and re-fitting the model iteratively.
#'
#' @param y A numeric matrix of observations, with each row corresponding to a time point.
#' @param init_par A numeric vector of initial parameter values for optimization.
#' @param build A function whose first argument is a parameter vector and that returns a \code{dlm} model (as used in \code{dlm::dlmMLE()}). The \code{specify_SSM} function can be used to create this \code{build} function. Additional named arguments can be declared in \code{build} and supplied via \code{build_args}.
#' @param build_args An optional named list of additional arguments to forward to \code{build} on every call. For example, \code{build_args = list(u = covariate_matrix)} allows \code{build} to be written as \code{function(par, u) \{ ... \}} without requiring a closure or factory function. Default is \code{list()}.
#' @param num_lambdas Integer. The number of \eqn{\lambda} values to evaluate. Ignored if \code{custom_lambdas} is specified. Default is 20.
#' @param custom_lambdas Optional numeric vector. If supplied, these are the exact \eqn{\lambda} values used for model fitting. If not provided or set to \code{NA}, then \code{num_lambdas} \eqn{\lambda}'s are automatically chosen.
#' @param cores Integer. Number of CPU cores to use for parallel processing. Default is 1 (sequential execution).
#' @param B Integer. Maximum number of iterations per \eqn{\lambda}. Default is 50.
#' @param lower Optional numeric vector of lower bounds for optimization. If \code{NA}, defaults to \code{-Inf} for all parameters. Must be of same length as \code{init_par}.
#' @param upper Optional numeric vector of upper bounds for optimization. If \code{NA}, defaults to \code{Inf} for all parameters. Must be of same length as \code{init_par}.
#' @param tol Tolerance level for checking convergence of the ROAMS procedure. Default is \code{1e-4}.
#' @param lambda_min Minimum \eqn{\lambda} value to consider when constructing the sequence of \eqn{\lambda}'s. Ignored if \code{custom_lambdas} is specified. Default is 2.
#' @param excessive_outliers_iter_limit Integer. Maximum number of iterations allowed where \eqn{\ge 50\%} of timepoints are flagged as outliers. This many outliers suggests \eqn{\lambda} is too low. Allows ROAMS to get through these \eqn{\lambda}'s quicker. Default is 1.
#' @param control A named list of control options to pass to \code{optim} via \code{dlm::dlmMLE()}. Default is \code{list(parscale = init_par)}, which can help the optimizer if parameters are on vastly different scales.
#' @param thresholding Character string specifying the outlier thresholding rule. Either \code{"hard"} (default) or \code{"soft"}. Hard thresholding fully removes flagged observations (sets them to missing), corresponding to an L0 penalty on the outlier shifts \eqn{\pmb{\gamma}}. Soft thresholding applies a continuous LASSO shrinkage to each \eqn{\pmb{\gamma}_t}, partially adjusting observations rather than discarding them entirely, which corresponds to an L1 penalty.
#'
#' @return If more than one \eqn{\lambda} values are used, returns an object of class \code{roams_SSM_list} — a list containing a \code{roams_SSM} model for each \eqn{\lambda}. If only one \eqn{\lambda} value is used (i.e. \code{custom_lambdas} is manually specified as a single value), returns a single \code{roams_SSM} object.
#'
#' Each \code{roams_SSM} object includes:
#' \itemize{
#'   \item \code{lambda} - The \eqn{\lambda} value used.
#'   \item \code{prop_outlying} - Proportion of non-missing time points identified as outliers.
#'   \item \code{BIC} - Bayesian Information Criterion of the final model.
#'   \item \code{loglik} - Log-likelihood of the fitted model.
#'   \item \code{RSS} - Residual sum of squares.
#'   \item \code{gamma} - Matrix of estimated outlier adjustments.
#'   \item \code{iterations} - Number of iterations performed.
#'   \item Optimization output from \code{dlm::dlmMLE()} from the final iteration.
#'   \item \code{y} - The original data matrix.
#'   \item \code{build} - The original build function used to specify the model.
#' }
#'
#' @details
#' The ROAMS procedure alternates between estimating model parameters via maximum likelihood and identifying outlying observations based on Mahalanobis distance of residuals. For each iteration:
#' \enumerate{
#'   \item A \code{dlm} model is fit using \code{dlm::dlmMLE()}.
#'   \item Mahalanobis distance of residuals (Mahalanobis residuals) are computed.
#'   \item Outlier shift estimates \eqn{\pmb{\gamma}_t} are updated using the selected thresholding rule.
#'   \item The adjusted data \code{adj_y} is updated and passed to the next iteration.
#' }
#'
#' Under \strong{hard thresholding}, observations whose penalised Mahalanobis score \eqn{d_t^2 + \log|\mathbf{S}_{t|t-1}|} exceeds \eqn{\lambda^2} are fully removed (set to missing), equivalent to setting \eqn{\hat{\pmb{\gamma}}_t = \mathbf{r}_t}. Under \strong{soft thresholding}, a LASSO shrinkage operator is applied: \eqn{\hat{\pmb{\gamma}}_t = \max(1 - \lambda / d_t,\, 0)\, \mathbf{r}_t}, where \eqn{d_t} is the Mahalanobis distance of the one-step-ahead residual. Partial adjustments are retained in the data (\code{adj_y = y - gamma}) rather than treating observations as missing.
#'
#' The algorithm stops when the change in parameters and outlier estimates is sufficiently small or if too many outliers are detected (more than 50\% of complete observations).
#'
#' @seealso \code{\link[dlm]{dlmMLE}}, \code{\link{best_BIC_model}}, \code{\link{outlier_target_model}}, \code{\link{get_attribute}}, \code{\link{autoplot.roams_SSM_list}}, \code{\link{attach_insample_info}}, \code{\link{oos_filter}}, \code{\link{specify_SSM}}
#'
#' @references She, Y., & Owen, A. B. (2011). Outlier Detection Using Nonconvex Penalized Regression. \emph{Journal of the American Statistical Association, 106}(494), 626–639. https://doi.org/10.1198/jasa.2011.tm10390
#'
#' @export
roams_SSM = function(
    y,
    init_par,
    build,
    build_args = list(),
    num_lambdas = 20,
    custom_lambdas = NA,
    cores = 1,
    B = 50,
    lower = NA,
    upper = NA,
    tol = 1e-4,
    lambda_min = 2,
    excessive_outliers_iter_limit = 1,
    control = list(parscale = init_par),
    thresholding = "hard"
    ) {

  if (ncol(y) > nrow(y)) {
    warning("Input data has more columns than rows. Did you forget to transpose your data? This function expects each row to represent a time point (i.e., observations in rows).")
  }

  if (!is.null(control$parscale) && any(control$parscale == 0)) {
    stop(
      "control$parscale contains zeros, which will cause optim() to fail. ",
      "This typically happens when init_par contains zeros, since the default ",
      "control sets parscale = init_par. ",
      "Either use non-zero initial values or supply a custom control, e.g. ",
      "control = list(parscale = rep(1, length(init_par)))."
    )
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
                         tol = tol,
                         excessive_outliers_iter_limit = excessive_outliers_iter_limit,
                         control = control,
                         thresholding = thresholding,
                         build_args = build_args)

    # Highest lambda is the supremum norm of mahalanobis residuals of classical fit.
    # classical$build is already wrapped, so pass it directly.
    highest_lambda = max(dlmInfo(y, y, classical, classical$build)$mahalanobis_residuals)
    lowest_lambda = lambda_min

    if (lowest_lambda >= highest_lambda) {
      stop("Automatic lambda selection failed. Lowest lambda = 1 is too large. Your data set may not have outliers. Try providing a custom lambda sequence via the 'custom_lambdas' argument.")
    }

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
                       tol = tol,
                       excessive_outliers_iter_limit = excessive_outliers_iter_limit,
                       control = control,
                       thresholding = thresholding,
                       build_args = build_args)

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
                           tol = tol,
                           excessive_outliers_iter_limit = excessive_outliers_iter_limit,
                           control = control,
                           thresholding = thresholding,
                           build_args = build_args)

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
    tol,
    excessive_outliers_iter_limit = excessive_outliers_iter_limit,
    control,
    thresholding = "hard",
    build_args = list()
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
                                 tol,
                                 excessive_outliers_iter_limit,
                                 control,
                                 thresholding,
                                 build_args)
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
                                                       tol,
                                                       excessive_outliers_iter_limit,
                                                       control,
                                                       thresholding,
                                                       build_args))

    future::plan(future::sequential)
  }

  class(model_list) = "roams_SSM_list"
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
    tol,
    excessive_outliers_iter_limit,
    control,
    thresholding = "hard",
    build_args = list()
) {

  build = wrap_build(build, build_args)

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
  prop_outlying = 0
  excessive_outliers_offences = 0

  for (j in 1:B) {
    if (j != 1) {theta_old = fit$par}
    if (prop_outlying >= 0.5) {
      excessive_outliers_offences = excessive_outliers_offences + 1
      if (excessive_outliers_offences >= excessive_outliers_iter_limit) {
        break
      }
    }

    fit = dlm::dlmMLE(
      adj_y,
      parm = par,
      build = build,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = control
      )

    # Update gammas
    info_output = dlmInfo(y, adj_y, fit, build)
    r = y - info_output$predicted_observations
    gamma_old = gamma
    gamma = matrix(0, nrow = n, ncol = dim_obs)

    if (thresholding == "hard") {
      # Hard thresholding: fully remove flagged observations (L0 penalty)
      objective_increase = info_output$mahalanobis_residuals^2 + log(info_output$det_S)
      # NA timepoints should have the lowest penalised score so they are never flagged
      objective_increase = ifelse(info_output$mahalanobis_residuals == 0,
                                  min(min(objective_increase), lambda) - 1,
                                  objective_increase)
      gamma[objective_increase > lambda^2,] = r[objective_increase > lambda^2,]
      which_nz = which(rowSums(abs(gamma)) != 0)

      # Too many outliers detected will cause instability
      if (length(which_nz) / n_complete >= 0.5) {
        # Re-update gamma to only keep best nz timepoints
        keep_threshold = n_complete / 2 + (n - n_complete)
        which_nz = which(rank(objective_increase) > keep_threshold)
        gamma = matrix(0, nrow = n, ncol = dim_obs)
        gamma[which_nz,] = r[which_nz,]
      }

      adj_y = y
      adj_y[which_nz,] = NA

    } else if (thresholding == "soft") {
      # Soft thresholding: LASSO shrinkage (L1 penalty)
      # gamma_t = max(1 - lambda / d_t, 0) * r_t, where d_t is Mahalanobis distance
      mah_resid = info_output$mahalanobis_residuals
      scale_factor = ifelse(mah_resid > 0, pmax(1 - lambda / mah_resid, 0), 0)
      gamma = r * scale_factor
      which_nz = which(scale_factor > 0)

      # Too many outliers detected will cause instability
      if (length(which_nz) / n_complete >= 0.5) {
        # Restrict to top 50% of timepoints by Mahalanobis residual
        keep_threshold = n_complete / 2 + (n - n_complete)
        keep_rows = which(rank(mah_resid) > keep_threshold)
        gamma = matrix(0, nrow = n, ncol = dim_obs)
        keep_scale = ifelse(mah_resid[keep_rows] > 0,
                            pmax(1 - lambda / mah_resid[keep_rows], 0), 0)
        gamma[keep_rows,] = r[keep_rows,] * keep_scale
        which_nz = keep_rows[keep_scale > 0]
      }

      # Partially adjust observations rather than treating them as missing
      adj_y = y - gamma

    } else {
      stop("Invalid 'thresholding' value. Must be \"hard\" or \"soft\".")
    }

    gap_gamma = max(abs(gamma - gamma_old))
    gap_theta = max(abs(fit$par - theta_old))

    nz = sum(rowSums(abs(gamma_old)) != 0)  # compute based on k-1 iter of gamma
    prop_outlying = nz / n_complete

    if ((gap_gamma < tol) && (gap_theta < tol)) {
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
    list(build = build),       # wrapped build (callable as build(par))
    list(build_args = build_args)
    )

  class(model) = "roams_SSM"
  return(model)
}

# dlm returns state vectors (not matrices) when dim_state = 1.
# This coerces them back to proper matrices so all downstream code
# can use consistent two-dimensional indexing.
ensure_state_matrix = function(x, dim_state) {
  if (!is.matrix(x)) matrix(x, ncol = dim_state) else x
}

# Wrap a build function to pre-supply named extra arguments.
# When build_args is empty, returns build unchanged — no closure overhead,
# and model$build remains the user's original function.
wrap_build = function(build, build_args) {
  force(build)  # evaluate promise now, before the closure captures it
  if (length(build_args) == 0) return(build)
  function(par, ...) do.call(build, c(list(par), build_args))
}

# Apply time-varying matrix elements from covariates (X/JFF or X/JGG) at time t.
# Returns the base matrix with indexed elements replaced by X[t, k].
apply_covariate_matrix = function(base_mat, J_mat, X, t) {
  if (is.null(J_mat) || is.null(X)) return(base_mat)
  out = base_mat
  nz = which(J_mat != 0, arr.ind = TRUE)
  for (k in seq_len(nrow(nz))) {
    out[nz[k, 1], nz[k, 2]] = X[t, J_mat[nz[k, 1], nz[k, 2]]]
  }
  out
}

# Compute list of S_t = A_t %*% P_{t|t-1} %*% t(A_t) + V for each t.
# When JFF is NULL, A is constant and a vectorised path is used.
compute_S_list = function(P_tt_1_list, A_base, JFF, X_cov, V) {
  if (is.null(JFF)) {
    purrr::map(P_tt_1_list, ~ A_base %*% . %*% t(A_base) + V)
  } else {
    purrr::map2(P_tt_1_list, seq_along(P_tt_1_list), ~ {
      A_t = apply_covariate_matrix(A_base, JFF, X_cov, .y)
      A_t %*% .x %*% t(A_t) + V
    })
  }
}

# Compute an n x dim_obs matrix of A_t %*% state_t for t = 1, ..., n.
# state_mat has n+1 rows (row 1 is t=0); rows 2:(n+1) correspond to t=1,...,n.
# When JFF is NULL a single matrix multiply is used; otherwise per-step loop.
compute_obs_from_states = function(state_mat, A_base, JFF, X_cov, n) {
  if (is.null(JFF)) {
    (state_mat[2:(n + 1), , drop = FALSE]) %*% t(A_base)
  } else {
    do.call(rbind, lapply(seq_len(n), function(t) {
      drop(apply_covariate_matrix(A_base, JFF, X_cov, t) %*% state_mat[t + 1, ])
    }))
  }
}

IPOD_oos_robust_filter = function(y, par, build, threshold, multiplier = 1) {

  SSM_specs = build(par)

  Phi_base = SSM_specs$GG
  Sigma_w  = SSM_specs$W
  A_base   = SSM_specs$FF
  Sigma_v  = SSM_specs$V
  x_tt     = SSM_specs$m0
  P_tt     = SSM_specs$C0

  # Covariate components (may be NULL for models without exogenous inputs)
  X_cov = SSM_specs$X
  JFF   = SSM_specs$JFF
  JGG   = SSM_specs$JGG

  n         = nrow(y)
  dim_obs   = ncol(y)
  dim_state = nrow(Phi_base)

  x_tt_1 = NA
  P_tt_1 = NA
  y_tt_1 = NA
  S_t    = NA

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
    # Resolve time-varying matrices at time t
    A   = apply_covariate_matrix(A_base,   JFF, X_cov, t)
    Phi = apply_covariate_matrix(Phi_base, JGG, X_cov, t)

    x_tt_1 = Phi %*% x_tt
    P_tt_1 = Phi %*% P_tt %*% t(Phi) + Sigma_w
    y_tt_1 = A %*% x_tt_1
    S_t    = A %*% P_tt_1 %*% t(A) + Sigma_v
    inv_S_t = solve(S_t)

    if (any(is.na(y[t,]))) {
      mahalanobis_residuals[t] = 0
      x_tt = x_tt_1
      P_tt = multiplier * P_tt_1
    } else {

      mahalanobis_residuals[t] = drop(sqrt(t(y[t,] - y_tt_1) %*% inv_S_t %*% (y[t,] - y_tt_1)))

      if (mahalanobis_residuals[t] <= threshold) {
        K_t  = P_tt_1 %*% t(A) %*% inv_S_t
        x_tt = x_tt_1 + K_t %*% (y[t,] - y_tt_1)
        P_tt = P_tt_1 - K_t %*% A %*% P_tt_1
      } else {
        x_tt = x_tt_1
        P_tt = multiplier * P_tt_1
        outliers_flagged[t] = 1
      }
    }

    filtered_states[t,]    = x_tt
    filtered_observations[t,] = A %*% x_tt
    predicted_states[t,]   = x_tt_1
    predicted_observations[t,] = y_tt_1
    filtered_states_var[[t]]   = P_tt
    predicted_states_var[[t]]  = P_tt_1
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
# Mainly used for getting residuals and predictions from ROAMS SSMs.
dlmInfo = function(y, adj_y, model, build) {

  filter_output   = dlm::dlmFilter(adj_y, mod = build(model$par))
  smoother_output = dlm::dlmSmooth(filter_output)

  A_base    = filter_output$mod$FF
  X_cov     = filter_output$mod$X    # NULL when no covariates
  JFF       = filter_output$mod$JFF  # NULL when no covariates
  V         = filter_output$mod$V
  n         = nrow(y)
  dim_state = nrow(filter_output$mod$GG)

  filter_output$m   = ensure_state_matrix(filter_output$m,   dim_state)
  smoother_output$s = ensure_state_matrix(smoother_output$s, dim_state)

  P_tt_1_list = dlm::dlmSvd2var(filter_output$U.R, filter_output$D.R)
  S     = compute_S_list(P_tt_1_list, A_base, JFF, X_cov, V)
  inv_S = purrr::map(S, .f = solve)
  mahalanobis_residuals = purrr::map2_dbl(
    apply(y - filter_output$f, 1, c, simplify = FALSE),
    inv_S,
    ~ drop(t(.x) %*% .y %*% .x)) |> sqrt()

  det_S = purrr::map_dbl(S, .f = det)
  mahalanobis_residuals = ifelse(is.na(mahalanobis_residuals), 0, mahalanobis_residuals)

  return(list(
    smoothed_observations  = compute_obs_from_states(smoother_output$s, A_base, JFF, X_cov, n),
    filtered_observations  = compute_obs_from_states(filter_output$m,   A_base, JFF, X_cov, n),
    predicted_observations = filter_output$f,
    mahalanobis_residuals  = mahalanobis_residuals,
    det_S = det_S
  ))

}

#' Attach In-Sample Information to Fitted State Space Model
#'
#' Attaches detailed in-sample information—such as predicted, filtered, and smoothed states and observations—to a model object fitted using any of the package’s supported SSM estimation methods.
#' These quantities are not stored by default in model objects due to their potentially large memory footprint.
#'
#' @param model A fitted model object of class \code{roams_SSM}, \code{classical_SSM}, \code{oracle_SSM}, \code{huber_robust_SSM}, or \code{trimmed_robust_SSM}.
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
#' For models of class \code{roams_SSM}, \code{classical_SSM}, or \code{oracle_SSM}, the following additional elements are also attached:
#' \describe{
#'   \item{\code{smoothed_states}}{Posterior means of hidden states using all data.}
#'   \item{\code{smoothed_observations}}{Posterior mean of the observed series based on smoothed states.}
#'   \item{\code{smoothed_states_var}}{List of smoothed state variance matrices.}
#' }
#'
#' These smoothed attributes are obtained using the RTS smoothing algorithm  (Rauch et al. 1965).
#'
#' @details
#' The attached outputs enable richer diagnostics, outlier inspection, and plotting.
#' For \code{huber_robust_SSM} and \code{trimmed_robust_SSM} models, in-sample information is computed using a custom robust filtering function, and smoothed quantities (\code{smoothed_states}, \code{smoothed_observations}, and \code{smoothed_states_var}) are \strong{not available}.
#' This function should only be applied once to a model object.
#'
#' @references Rauch, H.E., Tung, F., Striebel, C.T. (1965). Maximum likelihood estimates of linear dynamic systems. \emph{AIAA Journal 3}(8), 1445–1450. https://doi.org/10.2514/3.3166
#'
#' @seealso \code{\link{oos_filter}}
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

  } else if (inherits(model, "roams_SSM")){
    adj_y = y
    adj_y[which(rowSums(abs(model$gamma)) != 0),] = NA
  } else if (inherits(model, "classical_SSM")) {
    adj_y = y
  } else if (inherits(model, "oracle_SSM")) {
    adj_y = y
    adj_y[model$outlier_locs != 0,] = NA
  } else {
    stop("Invalid model class. Expected 'roams_SSM' or 'classical_SSM' or 'oracle_SSM' or 'huber_robust_SSM' or 'trimmed_robust_SSM'.")
  }

  filter_output   = dlm::dlmFilter(adj_y, mod = model$build(model$par))
  smoother_output = dlm::dlmSmooth(filter_output)

  A_base    = filter_output$mod$FF
  X_cov     = filter_output$mod$X    # NULL when no covariates
  JFF       = filter_output$mod$JFF  # NULL when no covariates
  V         = filter_output$mod$V
  n         = nrow(y)
  dim_state = nrow(filter_output$mod$GG)

  filter_output$m   = ensure_state_matrix(filter_output$m,   dim_state)
  filter_output$a   = ensure_state_matrix(filter_output$a,   dim_state)
  smoother_output$s = ensure_state_matrix(smoother_output$s, dim_state)

  P_tt_1      = dlm::dlmSvd2var(filter_output$U.R, filter_output$D.R)
  S           = compute_S_list(P_tt_1, A_base, JFF, X_cov, V)
  P_tt        = dlm::dlmSvd2var(filter_output$U.C, filter_output$D.C)
  P_tt_smooth = dlm::dlmSvd2var(smoother_output$U.S, smoother_output$D.S)

  inv_S = purrr::map(S, ~ solve(.))
  mahalanobis_residuals = purrr::map2_dbl(
    apply(y - filter_output$f, 1, c, simplify = FALSE),
    inv_S,
    ~ drop(t(.x) %*% .y %*% .x)) |> sqrt()

  mahalanobis_residuals = ifelse(is.na(mahalanobis_residuals), 0, mahalanobis_residuals)

  output = c(model,
    list(
    smoothed_states  = smoother_output$s[2:(n + 1), , drop = FALSE],
    filtered_states  = filter_output$m[2:(n + 1), , drop = FALSE],
    predicted_states = filter_output$a,
    smoothed_observations  = compute_obs_from_states(smoother_output$s, A_base, JFF, X_cov, n),
    filtered_observations  = compute_obs_from_states(filter_output$m,   A_base, JFF, X_cov, n),
    predicted_observations = filter_output$f,
    smoothed_states_var = P_tt_smooth,
    filtered_states_var = P_tt,
    predicted_states_var = P_tt_1,
    predicted_observations_var = S,
    mahalanobis_residuals = mahalanobis_residuals
  ))

  if (inherits(model, "roams_SSM")) {
    class(output) = c("roams_SSM", "insample_info")
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
#' @param model A fitted model object of class \code{roams_SSM}, \code{classical_SSM}, \code{oracle_SSM}, \code{huber_robust_SSM}, or \code{trimmed_robust_SSM}.
#' @param build A function whose first argument is a parameter vector and that returns a \code{dlm} model object. The \code{specify_SSM} function can be used to create this \code{build} function.
#' @param build_args An optional named list of additional arguments to forward to \code{build}. Should match what was supplied to the fitting function (e.g. \code{roams_SSM}) unless a different covariate matrix is needed for the out-of-sample period. Default is \code{list()}.
#' @param outlier_locs A logical or binary vector of the same length as \code{nrow(y)}, indicating time points to be treated as missing (i.e., time points that are known to be outliers). Used only with \code{oracle_SSM} models.
#' @param threshold Mahalanobis distance threshold for detecting out-of-sample outliers in \code{roams_SSM} models. Set to \code{Inf} to recover the usual Kalman filter. Default is \code{sqrt(qchisq(0.99, ncol(y)))}.
#' @param multiplier Multiplier for how quickly the filter grows its filtered state variance (uncertainty) after detecting an outlier in \code{roams_SSM} models. It is the tuneable parameter \eqn{b} of the `fast-updating threshold' filter. Only works if \code{threshold} is not \code{Inf}. Default is \code{2}.
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
#'   \item{\code{outliers_flagged}}{Vector of 1's and 0's indicating whether timepoints are flagged as outlying or not based on the \code{threshold} supplied (only available if \code{model} is of class \code{roams_SSM}).}
#' }
#'
#' @details
#' The function reuses the model's fitted parameters to generate inference on new data \code{y_oos}. Robust variants use appropriate robust filters, while the classical and oracle models use standard Kalman filtering. For \code{oracle_SSM} models, observations flagged in \code{outlier_locs} are treated as missing during filtering.
#'
#' @seealso \code{\link{attach_insample_info}}, \code{\link{specify_SSM}}
#'
#' @export
oos_filter = function(y_oos, model, build,
                      build_args = list(),
                      outlier_locs = rep(0, nrow(y_oos)),
                      threshold = sqrt(qchisq(0.99, ncol(y_oos))),
                      multiplier = 2) {

  y     = y_oos
  build = wrap_build(build, build_args)

  if (inherits(model, "huber_robust_SSM")) {
    oos_output = ruben_filter(model$par, y, build, obj_type = "huber")
    return(oos_output)
  } else if (inherits(model, "trimmed_robust_SSM")) {
    oos_output = ruben_filter(model$par, y, build, obj_type = "trimmed")
    return(oos_output)

  } else if (inherits(model, "roams_SSM")){

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
    stop("Invalid model class. Expected 'roams_SSM' or 'classical_SSM' or 'oracle_SSM' or 'huber_robust_SSM' or 'trimmed_robust_SSM'.")
  }

  filter_output = dlm::dlmFilter(adj_y, mod = build(model$par))

  A_base    = filter_output$mod$FF
  X_cov     = filter_output$mod$X    # NULL when no covariates
  JFF       = filter_output$mod$JFF  # NULL when no covariates
  V         = filter_output$mod$V
  n         = nrow(y)
  dim_state = nrow(filter_output$mod$GG)

  filter_output$m = ensure_state_matrix(filter_output$m, dim_state)
  filter_output$a = ensure_state_matrix(filter_output$a, dim_state)

  P_tt_1 = dlm::dlmSvd2var(filter_output$U.R, filter_output$D.R)
  S      = compute_S_list(P_tt_1, A_base, JFF, X_cov, V)
  P_tt   = dlm::dlmSvd2var(filter_output$U.C, filter_output$D.C)

  inv_S = purrr::map(S, ~ solve(.))
  mahalanobis_residuals = purrr::map2_dbl(
    apply(y - filter_output$f, 1, c, simplify = FALSE),
    inv_S,
    ~ drop(t(.x) %*% .y %*% .x)) |> sqrt()

  mahalanobis_residuals = ifelse(is.na(mahalanobis_residuals), 0, mahalanobis_residuals)

  return(list(
             filtered_states  = filter_output$m[2:(n + 1), , drop = FALSE],
             predicted_states = filter_output$a,
             filtered_observations  = compute_obs_from_states(filter_output$m, A_base, JFF, X_cov, n),
             predicted_observations = filter_output$f,
             filtered_states_var = P_tt,
             predicted_states_var = P_tt_1,
             predicted_observations_var = S,
             mahalanobis_residuals = mahalanobis_residuals
           ))

}



