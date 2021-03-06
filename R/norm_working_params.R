#' Transform normal HMM parameters from natural to working
#'
#' The function transforms the natural normal HMM parameters that have
#' additional constraints into working parameters that incorporate the
#' constraints.
#'
#' @param num_states The number of states in the desired HMM.
#' @param num_variables The number of variables in the data.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param mu A list of matrices containing the means of the state dependent
#'   normal distribution. Each matrix corresponds to a different variable,
#'   each row corresponds to a different subject and each column corresponds
#'   to a different state.
#' @param sigma A list of matrices containing the standard deviations of the
#'  state dependent normal distribution. Each matrix corresponds to a
#'  different variable, each row corresponds to a different subject and each
#'  column corresponds to a different state.
#' @param beta A matrix of regression coefficients for the effect of the
#'   covariates on the transition probability matrix `gamma`.
#' @param delta A list of the initial state distributions for each subject.
#'
#' @return A single vector containing working parameters.
#' @export
#' @examples
#' # define values of parameters
#' num_states    <- 2
#' num_variables <- 2
#' num_subjects  <- 2
#' mu            <- list(matrix(c(1, 5, 2, 4), 2, 2, byrow = TRUE),
#'                       matrix(c(1, 5, 2, 4), 2, 2, byrow = TRUE))
#' sigma         <- list(matrix(c(1, 2, 1, 1.5), 2, 2, byrow = TRUE),
#'                       matrix(c(1, 2, 1, 1.5), 2, 2, byrow = TRUE))
#' beta          <- matrix(c(-2, 0, 0), nrow = 1, ncol = 3)
#' delta         <- list(c(1/2, 1/2), c(1/2, 1/2))
#'
#' #transform to working parametes
#' norm_working_params(num_states, num_variables, num_subjects,
#'                     mu, sigma, beta, delta)

norm_working_params <- function(num_states, num_variables, num_subjects,
                                mu, sigma, beta, delta) {
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
