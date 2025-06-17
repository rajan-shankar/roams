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
best_BIC_model = function(model_list) {
  valid_indexes = which(get_attribute(model_list, "prop_outlying") < 0.5)
  model_list = model_list[valid_indexes]
  class(model_list) = "robularized_SSM_list"

  best_index = which.min(get_attribute(model_list, "BIC"))
  return(model_list[[best_index]])
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
outlier_target_model = function(model_list, target) {
  distances = abs(target - get_attribute(model_list, "prop_outlying"))
  index = which.min(distances)
  return(model_list[[index]])
}

