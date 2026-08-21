#' Select Best Model Based on an Information Criterion
#'
#' Extracts the best-fitting model from a \code{roams_SSM_list} object according to
#' a chosen information criterion, while excluding models with more than 50\% outlying observations.
#'
#' The criterion is computed as \eqn{-2\ell + k \cdot p}, where \eqn{\ell} is the
#' log-likelihood, \eqn{p} is the number of non-zero components, and \eqn{k} is
#' the penalty multiplier. The \code{ic} argument provides named shortcuts; \code{k}
#' allows a fully custom penalty and overrides \code{ic} when non-\code{NULL}.
#'
#' @param model_list An object of class \code{roams_SSM_list}.
#' @param ic Character string specifying the information criterion. One of
#'   \code{"BIC"} (default), \code{"AIC"}, or \code{"HQIC"} (Hannan-Quinn).
#'   Ignored if \code{k} is supplied.
#' @param k Optional numeric penalty multiplier per model parameter, as in
#'   \code{\link[stats]{step}}. Overrides \code{ic} when non-\code{NULL}.
#'   Common values: \code{2} (AIC), \code{log(n)} (BIC), \code{2 * log(log(n))} (HQIC).
#'
#' @return A single \code{roams_SSM} object corresponding to the model with the smallest
#'   value of the chosen information criterion among those with fewer than 50\% outlying observations.
#'
#' @examples
#' \dontrun{
#' # BIC (default)
#' best_BIC_model(models)
#'
#' # Named shortcuts
#' best_BIC_model(models, ic = "AIC")
#' best_BIC_model(models, ic = "HQIC")
#'
#' # Custom penalty via k (overrides ic)
#' best_BIC_model(models, k = 2)
#' }
#'
#' @export
best_BIC_model = function(model_list, ic = c("BIC", "AIC", "HQIC"), k = NULL) {
  ic = match.arg(ic)

  valid_indexes = which(get_attribute(model_list, "prop_outlying") < 0.5)
  model_list = model_list[valid_indexes]
  class(model_list) = "roams_SSM_list"

  loglik = get_attribute(model_list, "loglik")
  nz = get_attribute(model_list, "gamma") |>
    purrr::map_dbl(\(g) sum(rowSums(abs(g)) != 0))

  if (!is.null(k)) {
    # User-supplied penalty overrides ic
    ic_values = k * nz - 2 * loglik
  } else {
    n_complete = get_attribute(model_list, "y") |>
      purrr::map_dbl(\(y) sum(!is.na(y)))
    ic_values = switch(ic,
      AIC  = 2 * nz - 2 * loglik,
      BIC  = log(n_complete) * nz - 2 * loglik,
      HQIC = 2 * log(log(n_complete)) * nz - 2 * loglik
    )
  }

  best_index = which.min(ic_values)
  return(model_list[[best_index]])
}

#' Select Model Based on Target Outlier Proportion
#'
#' Extracts the model from a \code{roams_SSM_list} object whose estimated outlier proportion
#' is closest to a user-specified target.
#'
#' @param model_list An object of class \code{roams_SSM_list}.
#' @param target A numeric value between 0 and 1 indicating the desired proportion of outlying observations.
#'
#' @return A single \code{roams_SSM} object whose estimated outlier proportion is closest to \code{target}.
#'
#' @examples
#' \dontrun{
#' # Select the model with an outlier proportion closest to 10%
#' target_model <- outlier_target_model(models, target = 0.1)
#' }
#'
#' @export
outlier_target_model = function(model_list, target) {
  distances = abs(target - get_attribute(model_list, "prop_outlying"))
  index = which.min(distances)
  return(model_list[[index]])
}

