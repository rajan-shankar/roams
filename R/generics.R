#' @export
print.robularized_SSM_list = function(model_list) {
  print("testing list")
}

#' @export
print.robularized_SSM = function(model) {
  print("testing")
}

# I HAVE MADE THIS NOT A GENERIC FOR NOW. THERE WERE ISSUES WITH FURRR AND PLAN.
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
autoplot = function(model_list, attribute = "BIC") {

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

