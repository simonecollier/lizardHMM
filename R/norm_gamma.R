norm_gamma <- function(num_states, beta, design_vec) {
  #' Compute transition probability matrix
  #'
  #' This function computes the time dependent transition probability matrix
  #' from the corresponding vector `design_vec`  and the matrix of regression
  #' coefficients `beta`.
  #'
  #' @inheritParams norm_working_params
  #' @param design_vec A vector from the design matrix with indicates the values
  #'   of the covariates.
  #'
  #' @return A matrix.
  #' @export
  #' @examples
  #' norm_gamma(2, matrix(c(-2, 0, 0)), c(1, 0, 25))

  eta       <- matrix(0, nrow = num_states, ncol = num_states)
  gamma     <- matrix(0, nrow = num_states, ncol = num_states)
  betarow   <- 1
  for (i in 1:num_states) {
    for (j in 1:num_states) {
      if (i != j) {
        eta[i, j] <- sum(beta[betarow, ]*design_vec)
        betarow = betarow + 1
      }
    }
  }
  for (i in 1:num_states) {
    for (j in 1:num_states) {
      gamma[i, j] <- exp(eta[i, j])/sum(exp(eta[i, ]))
    }
  }
  gamma
}
