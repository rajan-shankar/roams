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
kalman_filter = function(y, par, build, gamma = array(0, dim(y))) {

  out = fn_filter(y = y,
                  par = par,
                  build = build,
                  gamma = gamma)
  return(out)
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
huber_robust_filter = function(y, par, build) {

  out = ruben_filter(y = y,
                     par = par,
                     build = build)

  return(out)
}

