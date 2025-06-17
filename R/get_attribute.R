#' Extract Attributes from a `robularized_SSM_list` or Single Model
#'
#' Retrieves a specified attribute from each model within a `robularized_SSM_list` object. If a single model of class `robularized_SSM` is provided, it is automatically wrapped into a list for extraction.
#'
#' @param model_list An object of class `robularized_SSM_list` or a single `robularized_SSM` model. May optionally include in-sample information added via [attach_insample_info()].
#' @param attribute A character string specifying the name of the attribute to extract. Must be one of the available scalar attributes (e.g. `"BIC"`, `"lambda"`) or list/matrix attributes (e.g. `"filtered_states"`, `"smoothed_states"`).
#'
#' @return
#' If `attribute` is a scalar attribute (e.g., `"BIC"` or `"lambda"`), returns a numeric or character vector containing that attribute across all models.
#'
#' If `attribute` is a list-like attribute (e.g., `"filtered_states"` or `"gamma"`), returns a list of those values, one per model.
#'
#' @details
#' Available attributes depend on whether in-sample information has been attached using [attach_insample_info()]. If not, only core model components (e.g. `par`, `gamma`, `y`) and scalar fitting metrics (e.g. `BIC`, `loglik`) are available.
#'
#' \strong{Vector Attributes (always available):}
#' \itemize{
#'   \item `"lambda"`: Regularization penalty value.
#'   \item `"prop_outlying"`: Proportion of outlying time points identified.
#'   \item `"BIC"`: Bayesian Information Criterion.
#'   \item `"loglik"`: Log-likelihood.
#'   \item `"RSS"`: Residual sum of squares.
#'   \item `"iterations"`: Number of optimization iterations.
#'   \item `"value"`: Final objective function value.
#'   \item `"counts"`: Optimization evaluation counts.
#'   \item `"convergence"`: Convergence status of optimizer.
#'   \item `"message"`: Optimizer termination message.
#' }
#'
#' \strong{List Attributes that are always available:}
#' \itemize{
#'   \item `"par"`, `"gamma"`, `"y"`
#' }
#'
#' \strong{List Attributes that are only available if in-sample info is attached:}
#' \itemize{
#'   \item `"smoothed_states"`, `"filtered_states"`, `"predicted_states"`
#'   \item `"smoothed_observations"`, `"filtered_observations"`, `"predicted_observations"`
#'   \item `"smoothed_states_var"`, `"filtered_states_var"`, `"predicted_states_var"`, `"predicted_observations_var"`
#'   \item `"mahalanobis_residuals"`
#' }
#'
#' If an invalid attribute is supplied, the function will list available attributes and suggest whether in-sample information might be missing.
#'
#' @seealso [attach_insample_info()], [robularized_SSM()]
#'
#' @examples
#' \dontrun{
#' models <- robularized_SSM(y, init_par, build)
#' lambdas <- get_attribute(models, "lambda")
#' BICs <- get_attribute(models, "BIC")
#'
#' # With in-sample info
#' models[[1]] <- attach_insample_info(models[[1]], build)
#' residuals <- get_attribute(models[[1]], "mahalanobis_residuals")
#' }
#'
#' @export
get_attribute = function(model_list, attribute) {

  if (inherits(model_list, "robularized_SSM")) {
    model_list = list(model_list)
    class(model_list) = "robularized_SSM_list"
  } else if (!inherits(model_list, "robularized_SSM_list")) {
    stop("'model_list' is not a 'robularized_SSM' or 'robularized_SSM_list'.")
  }

  if (inherits(model_list[[1]], "insample_info")) {
    list_attributes = c(
      "par",
      "gamma",
      "y",
      "smoothed_states",
      "filtered_states",
      "predicted_states",
      "smoothed_observations",
      "filtered_observations",
      "predicted_observations",
      "smoothed_states_var",
      "filtered_states_var",
      "predicted_states_var",
      "predicted_observations_var",
      "mahalanobis_residuals"
    )
  } else {
    list_attributes = c(
      "par",
      "gamma",
      "y"
    )
  }

  vector_attributes = c(
    "lambda",
    "prop_outlying",
    "BIC",
    "loglik",
    "RSS",
    "iterations",
    "value",
    "counts",
    "convergence",
    "message"
  )

  if (attribute %in% vector_attributes) {
    attributes = sapply(model_list, function(x) x[[attribute]])
    return(attributes)
  } else if (attribute %in% list_attributes) {
    attributes = lapply(model_list, function(x) x[[attribute]])
    return(attributes)
  } else {
    if (inherits(model_list[[1]], "insample_info")) {
      stop("This attribute does not exist. Available attributes are: \n",
           paste(c(vector_attributes, list_attributes), collapse = ", "), "\n")
    } else {
      stop("This attribute does not exist. Did you forget to attach_insample_info() to your model? Available attributes are: \n",
           paste(c(vector_attributes, list_attributes), collapse = ", "), "\n")
    }
  }

}
