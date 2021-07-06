norm_working_params <- function(num_states, num_variables, num_subjects,
                                mu, sigma, beta, delta) {
  #' Transform normal HMM parameters from natural to working
  #'
  #' The `norm_working_params` function transforms the natural normal HMM parameters that have
  #' additional constraints into working parameters that incorporate the
  #' constraints. The output is a single vector which includes all the working
  #' parameters.
  #'
  #' @param num_states The number of states in the HMM.
  #' @param num_subjects The number of subjects that generated the data being
  #'   fit with the HMM.
  #' @param num_variables The number of variables in the HMM data.
  #' @param mu The means for the normally distributed state dependent
  #'   distributions of the HMM. mu is a list of matrices with each matrix
  #'   corresponding to a different variable in the data being fit. The columns
  #'   of the matrices correspond to the state number and the rows correspond to
  #'   the subject number.
  #' @param sigma The standard deviations for the normally distributed state
  #'   dependent distributions of the HMM. sigma is a list of matrices with
  #'   each matrix corresponding to a different variable in the data being fit.
  #'   The columns of the matrices correspond to the state number and the rows
  #'   correspond to the subject number.
  #' @param beta The matrix of regression coefficients for the covariates
  #'   included in the estimated of the transition probability matrix.
  #' @param delta A list with each element being the initial state distribution
  #'   vector of the HMM for the subject corresponding to that index.

  tmu     <- numeric()
  tsigma  <- numeric()
  for (j in 1:num_variables) {
    tmu    <- c(tmu, as.vector(t(mu[[j]])))
    tsigma <- c(tsigma, log(as.vector(t(sigma[[j]]))))
  }
  if (num_states == 1) {
    return(tmu, tsigma)
  }
  tbeta  <- as.vector(beta)
  tdelta <- numeric()
  for (i in 1:num_subjects) {
    tdelta <- c(tdelta, log(delta[[i]][-1]/delta[[i]][1]))
  }
  c(tmu, tsigma, tbeta, tdelta)
}
