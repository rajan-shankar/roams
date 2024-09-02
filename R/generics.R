#' @export
print.robularized_SSM_list = function(model_list) {
  print("testing list")
}

#' @export
print.robularized_SSM = function(model) {
  print("testing")
}

#' @export
autoplot.robularized_SSM_list = function(model_list, attribute = "BIC") {

  vector_attributes = c(
    "lambda",
    "BIC",
    "RSS",
    "prop_outlying",
    "iterations"
  )

  if (!(attribute %in% vector_attributes)) {
    stop("This attribute does not exist.")
  }

  data = data.frame(
    lambda = get_attribute(model_list, "lambda"),
    BIC = get_attribute(model_list, "BIC"),
    attribute = get_attribute(model_list, attribute))

  data %>%
    ggplot() +
    aes(x = lambda, y = attribute) +
    geom_line(linewidth = 1) +
    geom_vline(data = . %>% dplyr::slice(which.min(BIC)),
               aes(xintercept = lambda),
               colour = "red", linetype = "dashed") +
    labs(x = latex2exp::TeX("$\\lambda$"),
         y = attribute) +
    theme_bw(base_size = 16)
}

