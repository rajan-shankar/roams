#' @export
print.robularized_SSM_list = function(model_list) {
  cat("Robularized SSM List with ", length(model_list), "models\n",
      "Use best_BIC_model() or outlier_target_model() to extract a single model.\n",
      "Use autoplot() to visualize the models.\n",
      "Use get_attribute() to extract attributes from the models.\n")
}

#' @export
print.robularized_SSM = function(model) {
  cat("Robularized SSM Model\n",
      "Lambda: ", model$lambda, "\n",
      "Outliers Detected: ", round(model$prop_outlying*100, 2), "%\n",
      "BIC: ", round(model$BIC, 3), "\n",
      "Log-Likelihood: ", round(model$loglik, 3), "\n",
      "RSS: ", round(model$RSS, 3), "\n",
      "IPOD Iterations: ", model$iterations, "\n",
      "Use $ to see more attributes.\n",)
}

#' @export
print.classical_SSM = function(model) {
  cat("Classical SSM Model\n",
      "Log-Likelihood: ", round(model$value, 3), "\n",
      "optim() Iterations: ", model$iterations, "\n",
      "Use $ to see more attributes.\n",)
}

#' @export
print.oracle_SSM = function(model) {
  cat("Oracle SSM Model\n",
      "Log-Likelihood: ", round(model$value, 3), "\n",
      "optim() Iterations: ", model$iterations, "\n",
      "Outlier Locations: ", model$outlier_locs, "n",
      "Use $ to see more attributes.\n",)
}

#' @export
print.huber_robust_SSM = function(model) {
  cat("Huber SSM Model\n",
      "Log-Likelihood: ", round(model$value, 3), "\n",
      "optim() Iterations: ", model$iterations, "\n",
      "Use $ to see more attributes.\n",)
}

#' @export
print.trimmed_robust_SSM = function(model) {
  cat("Trimmed SSM Model\n",
      "Log-Likelihood: ", round(model$value, 3), "\n",
      "optim() Iterations: ", model$iterations, "\n",
      "Alpha: ", model$alpha, "\n",
      "Use $ to see more attributes.\n",)
}

#' Autoplot for Robularized State Space Model Grid
#'
#' Generates a diagnostic plot for a list of robust state space models fit across a grid of \eqn{\lambda} values. The plot displays the specified model attribute (e.g., BIC, log-likelihood, proportion of outliers) against \eqn{\lambda}, with an optional vertical dashed line indicating the model with the lowest BIC among those with fewer than 45% outliers.
#'
#' @param model_list An object of class `"robularized_SSM_list"` as returned by [robularized_SSM()] when multiple \eqn{\lambda} values are used.
#' @param attribute A character string indicating which model attribute to plot. Options include `"lambda"`, `"prop_outlying"`, `"BIC"`, `"loglik"`, `"RSS"`, `"iterations"`, `"value"`, and `"counts"`. Defaults to `"BIC"`.
#'
#' @return A `ggplot` object showing the trajectory of the specified attribute across \eqn{\lambda}.
#'
#' @details
#' The red dashed vertical line indicates the model with the lowest BIC among models with less than 45% outlying time points, as a heuristic for robust model selection.
#'
#' @seealso [robularized_SSM()], [get_attribute()]
#'
#' @import ggplot2
#' @importFrom dplyr filter slice
#' @importFrom latex2exp TeX
#' @export
#' @method autoplot robularized_SSM_list
autoplot.robularized_SSM_list = function(model_list, attribute = "BIC") {

  vector_attributes = c(
    "lambda",
    "prop_outlying",
    "BIC",
    "loglik",
    "RSS",
    "iterations",
    "value",
    "counts"
  )

  if (!(attribute %in% vector_attributes)) {
    stop("This attribute does not exist or is not numeric.")
  }

  data = data.frame(
    lambda = get_attribute(model_list, "lambda"),
    BIC = get_attribute(model_list, "BIC"),
    prop_outlying = get_attribute(model_list, "prop_outlying"),
    attribute = get_attribute(model_list, attribute))

  data %>%
    ggplot() +
    aes(x = lambda, y = attribute) +
    geom_line(linewidth = 1) +
    geom_vline(data = . %>%
                 dplyr::filter(prop_outlying < 0.45) %>%
                 dplyr::slice(which.min(BIC)),
               aes(xintercept = lambda),
               colour = "red", linetype = "dashed") +
    labs(x = latex2exp::TeX("$\\lambda$"),
         y = attribute) +
    theme_bw(base_size = 16)
}

