#' Compute transition probability matrices
#'
#' This function computes the time dependent transition probability matrix
#' from the corresponding vector `design_vec` and the matrix of regression
#' coefficients `beta`.
#'
#' @param num_states The number of states in the desired HMM.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_time A value indicating the length of the data.
#' @param beta A matrix of regression coefficients for the effect of the
#'   covariates on the transition probability matrices `gamma`.
#' @param design A list of design matrices, one for each subject which
#' indicate the values of the each of the covariates (column) at each point
#' in time (row).
#'
#' @return A list containing a list for each subject/trial of the transition
#'   probability matrices at each time step.
#' @export
#' @example
#' # define parameters
#' num_states   <- 3
#' num_subjects <- 1
#' num_time     <- 1000
#' beta             <- matrix(c(0.01, 0.02, 0.001,
#'                              0.01, 0.03, 0.004,
#'                              0.01, 0.01, 0.003,
#'                              0.01, 0.04, 0.002,
#'                              0.01, 0.01, 0.004,
#'                              0.01, 0.03, 0.001),
#'                              ncol = 3, nrow = 6, byrow = TRUE)
#' design           <- list(matrix(0, nrow = 1000, ncol = 3))
#' design[[1]][, 1] <- 1 # First column is the intercept
#' design[[1]][, 2] <- sample(c(0, 1), size = 1000,
#'                            prob = c(0.3, 0.7), replace = TRUE)
#' design[[1]][, 3] <- rnorm(1000, mean = 5, sd = 1)
#'
#' # compute the transition probability matrices
#' fit_tpm(num_states, num_subjects, num_time, beta, design)
fit_tpm <- function(num_states, num_subjects, num_time, beta, design) {
  gamma <- list()
  for (i in 1:num_subjects) {
    gamma[[i]] <- list()
    for (t in 1:num_time) {
      eta             <- matrix(0, nrow = num_states, ncol = num_states)
      gamma[[i]][[t]] <- matrix(0, nrow = num_states, ncol = num_states)
      betarow         <- 1
      for (k in 1:num_states) {
        for (l in 1:num_states) {
          if (k != l) {
            eta[k, l] <- sum(beta[betarow, ]*design[[i]][t, ])
            betarow   <- betarow + 1
          }
        }
      }
      for (k in 1:num_states) {
        eta_row <- sum(exp(eta[k, ]))
        for (l in 1:num_states) {
          gamma[[i]][[t]][k, l] <- exp(eta[k, l])/eta_row
        }
      }
    }
  }
  gamma
}
