#' @export
print.roams_SSM_list = function(x, ...) {
  model_list = x
  cat("ROAMS SSM List with ", length(model_list), " models\n",
      " * Use best_BIC_model() or outlier_target_model() to extract a preferred model.\n",
      "   - Or, simply use list indexing such as [[1]] to extract a single model.\n",
      " * Use autoplot() to visualize the models.\n",
      " * Use get_attribute() to extract attributes from the models.\n", sep = "")
}

#' @export
print.roams_SSM = function(x, ...) {
  model = x
  cat("ROAMS SSM Model\n",
      " * Lambda: ", round(model$lambda, 3), "\n",
      " * Outliers Detected: ", round(model$prop_outlying*100, 2), "%\n",
      " * BIC: ", round(model$BIC, 3), "\n",
      " * Log-Likelihood: ", round(model$loglik, 3), "\n",
      " * Iterations: ", model$iterations, "\n",
      ifelse(inherits(model, "insample_info"),
             "In-sample information attached. Use $ to see more attributes.\n",
             "Use $ to see more attributes.\n"),
      sep = "")
}

#' @export
print.classical_SSM = function(x, ...) {
  model = x
  cat("Classical SSM Model\n",
      " * Log-Likelihood: ", round(model$value, 3), "\n",
      " * optim() Iterations: ", model$counts[1], "\n",
      ifelse(inherits(model, "insample_info"),
             "In-sample information attached. Use $ to see more attributes.\n",
             "Use $ to see more attributes.\n"),
      sep = "")
}

#' @export
print.oracle_SSM = function(x, ...) {
  model = x
  cat("Oracle SSM Model\n",
      " * Log-Likelihood: ", round(model$value, 3), "\n",
      " * optim() Iterations: ", model$counts[1], "\n",
      " * True Outlier Locations: ", model$outlier_locs, "\n",
      ifelse(inherits(model, "insample_info"),
             "In-sample information attached. Use $ to see more attributes.\n",
             "Use $ to see more attributes.\n"),
      sep = "")
}

#' @export
print.huber_robust_SSM = function(x, ...) {
  model = x
  cat("Huber Robust SSM Model\n",
      " * Log-Likelihood: ", round(model$value, 3), "\n",
      " * optim() Iterations: ", model$counts[1], "\n",
      ifelse(inherits(model, "insample_info"),
             "In-sample information attached. Use $ to see more attributes.\n",
             "Use $ to see more attributes.\n"),
      sep = "")
}

#' @export
print.trimmed_robust_SSM = function(x, ...) {
  model = x
  cat("Trimmed Robust SSM Model\n",
      " * Log-Likelihood: ", round(model$value, 3), "\n",
      " * optim() Iterations: ", model$counts[1], "\n",
      " * Alpha: ", model$alpha, "\n",
      ifelse(inherits(model, "insample_info"),
             "In-sample information attached. Use $ to see more attributes.\n",
             "Use $ to see more attributes.\n"),
      sep = "")
}

#' Autoplot for a ROAMS SSM List
#'
#' Generates diagnostic plots for a list of robust state space models fit across a sequence of \eqn{\lambda} values. Two model attributes are plotted against \eqn{\lambda}, each in its own panel, with a vertical dashed line indicating the model with the lowest BIC among those with fewer than 50\% outliers.
#'
#' @param object An object of class \code{roams_SSM_list} as returned by \code{\link{roams_SSM}} when multiple \eqn{\lambda} values are used.
#' @param attribute1 A character string indicating the first model attribute to plot on the top panel. Options include \code{"lambda"}, \code{"prop_outlying"}, \code{"BIC"}, \code{"loglik"}, \code{"RSS"}, \code{"iterations"}, \code{"value"}, and \code{"counts"}. Defaults to \code{"BIC"}.
#' @param attribute2 A character string indicating the second model attribute to plot on the bottom panel. Uses the same options as \code{attribute1}. Defaults to \code{"prop_outlying"}.
#' @param ... Other arguments passed to specific methods. Not used in this method.
#'
#' @return A \code{ggplot} object arranged using the \pkg{patchwork} package, showing the specified attributes plotted against \eqn{\lambda} in two vertically stacked panels.
#'
#' @details
#' In each panel, a red dashed vertical line indicates the model with the lowest BIC among those with fewer than 50\% outlying time points, serving as a heuristic for robust model selection.
#'
#' @seealso \code{\link{roams_SSM}}, \code{\link{get_attribute}}
#'
#' @importFrom ggplot2 autoplot
#'
#' @export
#' @method autoplot roams_SSM_list

autoplot.roams_SSM_list = function(object,
                                   attribute1 = "BIC",
                                   attribute2 = "prop_outlying",
                                   ...) {

  model_list = object

  vector_attributes = c(
    "lambda",
    "prop_outlying",
    "BIC",
    "loglik",
    "RSS",
    "iterations",
    "value"
  )

  if (!(attribute1 %in% vector_attributes) |
      !(attribute2 %in% vector_attributes)) {
    stop("At least one of the attributes specified does not exist or is not numeric.")
  }

  data = data.frame(
    lambda = get_attribute(model_list, "lambda"),
    BIC = get_attribute(model_list, "BIC"),
    prop_outlying = get_attribute(model_list, "prop_outlying"),
    attribute1 = get_attribute(model_list, attribute1),
    attribute2 = get_attribute(model_list, attribute2)
    )

  # Compute line position (same for both panels)
  vline_data = dplyr::filter(data, prop_outlying < 0.5) |>
    dplyr::slice(which.min(BIC))

  p1 = ggplot2::ggplot(data, ggplot2::aes(x = lambda, y = attribute1)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_vline(
      data = vline_data,
      ggplot2::aes(xintercept = lambda),
      colour = "red", linetype = "dashed"
    ) +
    ggplot2::labs(
      x = latex2exp::TeX("$\\lambda$"),
      y = attribute1
    ) +
    ggplot2::theme_bw(base_size = 16)

  if (attribute1 == "prop_outlying") {
    p1 = p1 + ggplot2::coord_cartesian(ylim = c(0, NA))
  }

  p2 = ggplot2::ggplot(data, ggplot2::aes(x = lambda, y = attribute2)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_vline(
      data = vline_data,
      ggplot2::aes(xintercept = lambda),
      colour = "red", linetype = "dashed"
    ) +
    ggplot2::labs(
      x = latex2exp::TeX("$\\lambda$"),
      y = attribute2
    ) +
    ggplot2::theme_bw(base_size = 16)

  if (attribute2 == "prop_outlying") {
    p2 = p2 + ggplot2::coord_cartesian(ylim = c(0, NA))
  }

  p1 + p2 + patchwork::plot_layout(
    nrow = 2, axes = "collect_x", axis_titles = "collect_x"
  )
}

