norm_gamma <- function(num_states, num_subjects, num_time, beta, design) {
  #' Compute transition probability matrix
  #'
  #' This function computes the time dependent transition probability matrix
  #' from the corresponding vector `design_vec`  and the matrix of regression
  #' coefficients `beta`.
  #'
  #' @inheritParams norm_working_params
  #' @param num_time A value indicating the length of the data.
  #' @param design A list of design matrices, one for each subject which
  #' indicate the values of the each of the covariates (column) at each point
  #' in time (row).
  #'
  #' @return A matrix.
  #' @export
  #' @examples
  #' norm_gamma(2, matrix(c(-2, 0, 0)), c(1, 0, 25))

  gamma <- list()
  for (i in num_subjects) {
    gamma[[i]] <- list()
    for (t in num_time) {
      eta             <- matrix(0, nrow = num_states, ncol = num_states)
      gamma[[i]][[t]] <- matrix(0, nrow = num_states, ncol = num_states)
      betarow         <- 1
      for (k in num_states) {
        for (l in num_states) {
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
