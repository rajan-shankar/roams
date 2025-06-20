#' Extract Attributes from a Robularized SSM List
#'
#' Retrieves a specified attribute from each model within a \code{robularized_SSM_list} object. Also works if a single model of class \code{robularized_SSM} is provided.
#'
#' @param model_list An object of class \code{robularized_SSM_list} or a single \code{robularized_SSM} model. May optionally include in-sample information added via \code{attach_insample_info}.
#' @param attribute A character string specifying the name of the attribute to extract. Must be one of the available scalar attributes (e.g. \code{BIC}, \code{lambda}) or list/vector attributes (e.g. \code{filtered_states}, \code{smoothed_states}).
#'
#' @return
#' If \code{attribute} is a scalar attribute (e.g., \code{BIC} or \code{lambda}), returns a numeric or character vector containing that attribute across all models.
#'
#' If \code{attribute} is a list-like attribute (e.g., \code{par} or \code{gamma}), returns a list of those values, across all models.
#'
#' @details
#' Available attributes depend on whether in-sample information has been attached using \code{attach_insample_info()}. If not, only core model components (e.g. \code{par}, \code{gamma}, \code{y}) and scalar metrics (e.g. \code{BIC}, \code{loglik}) are available.
#'
#' \strong{Scalar Attributes (always available):}
#' \itemize{
#'   \item \code{lambda}: Regularization penalty value.
#'   \item \code{prop_outlying}: Proportion of outlying time points identified.
#'   \item \code{BIC}: Bayesian Information Criterion.
#'   \item \code{loglik}: Log-likelihood.
#'   \item \code{RSS}: Residual sum of squares.
#'   \item \code{iterations}: Number of IPOD iterations.
#'   \item \code{value}: Final objective function value.
#'   \item \code{counts}: Optimization evaluation counts.
#'   \item \code{convergence}: Convergence status of optimizer.
#'   \item \code{message}: Optimizer termination message.
#' }
#'
#' \strong{List/Vector Attributes that are always available:}
#' \itemize{
#'   \item \code{par}, \code{gamma}, \code{y}
#' }
#'
#' \strong{List/Vector Attributes that are only available if in-sample info is attached:}
#' \itemize{
#'   \item \code{smoothed_states}, \code{filtered_states}, \code{predicted_states}
#'   \item \code{smoothed_observations}, \code{filtered_observations}, \code{predicted_observations}
#'   \item \code{smoothed_states_var}, \code{filtered_states_var}, \code{predicted_states_var}, \code{predicted_observations_var}
#'   \item \code{mahalanobis_residuals}
#' }
#'
#' Note that these 'in-sample info' attributes are typically only available if \code{model_list} is a single \code{robularized_SSM} with in-sample information already attached using \code{attach_insample_info()}. If \code{model_list} is a \code{robularized_SSM_list}, these attributes will not be available unless all models in the list have in-sample information attached.
#'
#' @seealso \code{\link{attach_insample_info}}, \code{\link{robularized_SSM}}
#'
#' @import tidyverse
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
