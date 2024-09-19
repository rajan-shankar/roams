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
compute_objective = function(y, par, build, gamma = array(0, dim(y))) {

  out = fn_filter(y = y,
                  par = par,
                  build = build,
                  gamma = gamma,
                  return_obj = TRUE)
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
compute_objective_no_gamma = function(y, par, build, outlier_locs) {

  out = no_gamma_oracle_filter(y = y,
                  par = par,
                  build = build,
                  outlier_locs = outlier_locs,
                  return_obj = TRUE)
  return(out)
}


