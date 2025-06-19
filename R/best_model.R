#' Select Best Model Based on BIC
#'
#' Extracts the best-fitting model from a \code{robularized_SSM_list} object according to
#' the Bayesian Information Criterion (BIC), while excluding models with more than 50\% outlying observations.
#'
#' @param model_list An object of class \code{robularized_SSM_list}
#'
#' @return A single \code{robularized_SSM} object corresponding to the model with the smallest BIC
#' among those with fewer than 50\% outlying observations.
#'
#' @examples
#' \dontrun{
#' # Assuming `models` is a robularized_SSM_list:
#' best_model <- best_BIC_model(models)
#' }
#'
#'
#' @export
best_BIC_model = function(model_list) {
  valid_indexes = which(get_attribute(model_list, "prop_outlying") < 0.5)
  model_list = model_list[valid_indexes]
  class(model_list) = "robularized_SSM_list"

  best_index = which.min(get_attribute(model_list, "BIC"))
  return(model_list[[best_index]])
}

#' Select Model Based on Target Outlier Proportion
#'
#' Extracts the model from a \code{robularized_SSM_list} object whose estimated outlier proportion
#' is closest to a user-specified target.
#'
#' @param model_list An object of class \code{robularized_SSM_list}.
#' @param target A numeric value between 0 and 1 indicating the desired proportion of outlying observations.
#'
#' @return A single \code{robularized_SSM} object whose estimated outlier proportion is closest to \code{target}.
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

