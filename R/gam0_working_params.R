#' Transform gamma HMM parameters from natural to working
#'
#' The function transforms the natural gamma HMM parameters that have
#' additional constraints into working parameters that incorporate the
#' constraints. Zero inflated gamma funciton.
#'
#' @param num_states The number of states in the desired HMM.
#' @param num_variables The number of variables in the data.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param alpha A list of matrices containing the shape parameters of the state
#'   dependent gamma distributions. Each matrix corresponds to a different
#'   variable, each row corresponds to a different subject and each column
#'   corresponds to a different state.
#' @param theta A list of matrices containing the scale parameters of the
#'  state dependent gamma distributions. Each matrix corresponds to a
#'  different variable, each row corresponds to a different subject and each
#'  column corresponds to a different state.
#' @param zweight A list of vectors, one for each variable, containing the
#'  weights for the zero inflated gamma distribution.
#' @param beta A matrix of regression coefficients for the effect of the
#'   covariates on the transition probability matrix `gamma`.
#' @param delta A list of the initial state distributions for each subject.
#'
#' @return A single vector containing working parameters.
#' @export
gam0_working_params <- function(num_states, num_variables, num_subjects,
                                alpha, theta, zweight, beta, delta) {
  talpha   <- numeric()
  ttheta   <- numeric()
  for (j in 1:num_variables) {
    talpha <- c(talpha,  log(as.vector(t(alpha[[j]]))))
    ttheta <- c(ttheta, log(as.vector(t(theta[[j]]))))
  }
  if (num_states == 1) {
    return(talpha, ttheta)
  }
  tzweight <- atanh(unlist(zweight))
  tbeta    <- as.vector(beta)
  tdelta   <- numeric()
  for (i in 1:num_subjects) {
    tdelta <- c(tdelta, log(delta[[i]][-1]/delta[[i]][1]))
  }
  c(talpha, ttheta, tzweight, tbeta, tdelta)
}

