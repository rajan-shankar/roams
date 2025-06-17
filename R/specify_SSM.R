#' A Cat Function
#'
#' This function allows you to express your love of cats.
#'
#' @param love Do you love cats? Defaults to TRUE.
#' @details
#' Additional details...
#' @returns description
#' @export
specify_SSM = function(
    state_transition_matrix,
    state_noise_var,
    observation_matrix,
    observation_noise_var,
    init_state_mean,
    init_state_var
) {

  return(list(
    "GG" = state_transition_matrix,
    "W" = state_noise_var,
    "FF" = observation_matrix,
    "V" = observation_noise_var,
    "m0" = init_state_mean,
    "C0" = init_state_var
  ))
}
