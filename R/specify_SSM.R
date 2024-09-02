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
    "state_transition_matrix" = state_transition_matrix,
    "state_noise_var" = state_noise_var,
    "observation_matrix" = observation_matrix,
    "observation_noise_var" = observation_noise_var,
    "init_state_mean" = init_state_mean,
    "init_state_var" = init_state_var
  ))
}
