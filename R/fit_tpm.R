#' Compute transition probability matrices
#'
#' This function computes the time dependent transition probability matrix
#' from the corresponding vector `design` and the matrix of regression
#' coefficients `beta`.
#'
#' @param num_states The number of states in the desired HMM.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_covariates The number of covariates in the data that the
#'   transition probability matrix depends on.
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
fit_tpm <- function(num_states, num_subjects, num_covariates, num_time,
                    beta, design) {
  gamma <- list()
  if (num_covariates == 0) {
    num_time <- 1
  }
  for (i in 1:num_subjects) {
    gamma[[i]] <- array(dim = c(num_states, num_states, num_time))
    for (t in 1:num_time) {
      eta               <- matrix(0, nrow = num_states, ncol = num_states)
      gamma[[i]][, , t] <- matrix(0, nrow = num_states, ncol = num_states)
      betarow           <- 1
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
          gamma[[i]][[k, l, t]] <- exp(eta[k, l])/eta_row
        }
      }
    }
  }
  if (num_covariates == 0) {
    for (i in num_subjects) {
      gamma[[i]] <- as.matrix(gamma[[i]][, , 1])
    }
  }
  gamma
}
