#' Specify a State-Space Model in DLM Format
#'
#' A helper function for users to construct a state-space model in the format expected by the \code{dlm} package, which is the package used internally for most of the model fitting functions.
#' This function returns a named list of model components which can be directly used in a user-defined
#' \code{build} function passed to any modeling function in this package.
#'
#' @param state_transition_matrix A square matrix specifying the state transition dynamics (\code{GG}).
#' @param state_noise_var A square matrix specifying the variance of the state noise (\code{W}).
#' @param observation_matrix A matrix mapping the state to the observations (\code{FF}).
#' @param observation_noise_var A square matrix specifying the variance of the observation noise (\code{V}).
#' @param init_state_mean A vector specifying the initial mean of the state (\code{m0}).
#' @param init_state_var A square matrix specifying the initial state covariance (\code{C0}).
#'
#' @details
#' The letters in the parentheses in the Arguments section correspond to the naming convention used in the \code{dlm} package.
#'
#' @return A named list with elements \code{GG}, \code{W}, \code{FF}, \code{V}, \code{m0}, and \code{C0},
#' suitable for use in a custom \code{build} function for modeling or for online filtering (e.g., using \code{\link{oos_filter}}).
#'
#' @examples
#' build_function = function(par) {
#'   phi_coef = par[1]
#'   Phi = diag(c(1+phi_coef, 1+phi_coef, 0, 0))
#'   Phi[1,3] = -phi_coef
#'   Phi[2,4] = -phi_coef
#'   Phi[3,1] = 1
#'   Phi[4,2] = 1
#'
#'   A = diag(4)[1:2,]
#'   Sigma_W = diag(c(par[2], par[3], 0, 0))
#'   Sigma_V = diag(c(par[4], par[5]))
#'
#'   mu0 = rep(0, 4)
#'   P0 = diag(rep(0, 4))
#'
#'   specify_SSM(
#'     state_transition_matrix = Phi,
#'     state_noise_var = Sigma_W,
#'     observation_matrix = A,
#'     observation_noise_var = Sigma_V,
#'     init_state_mean = mu0,
#'     init_state_var = P0
#'   )
#' }
#'
#' @seealso \code{\link{robularized_SSM}}, \code{\link{oos_filter}}
#'
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
