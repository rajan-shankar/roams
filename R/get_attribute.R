#' A Cat Function
#'
#' This function allows you to express your love of cats.
#'
#' @param love Do you love cats? Defaults to TRUE.
#' @details
#' Additional details...
#' @returns description
#' @export
get_attribute = function(model_list, attribute) {

  if (class(model_list) == "robularized_SSM") {
    model_list = list(model_list)
    class(model_list) = "robularized_SSM_list"
  } else if (class(model_list) != "robularized_SSM_list") {
    stop("model_list is not an output from the robularized_SSM() function.")
  }

  vector_attributes = c(
    "lambda",
    "BIC",
    "RSS",
    "prop_outlying",
    "iterations"
  )

  list_attributes = c(
    "par",
    "gamma",
    "filtered_states",
    "filtered_observations",
    "predicted_states",
    "predicted_observations",
    "predicted_observations_var",
    "mahalanobis_residuals"
  )

  if (attribute %in% vector_attributes) {
    attributes = sapply(model_list, function(x) x[[attribute]])
    return(attributes)
  } else if (attribute %in% list_attributes) {
    attributes = lapply(model_list, function(x) x[[attribute]])
    return(attributes)
  } else {
    stop("This attribute does not exist.")
  }

}
