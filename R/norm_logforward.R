#' Compute log forward probabilities
#'
#' This function computes the log forward probabilities of the data based on
#' the HMM hmm.
#'
#' @param x The data to be fit with an HMM in the form of a 3D array. The
#'   first index (row) corresponds to time, the second (column) to the
#'   variable number, and the third (matrix number) to the subject number.
#' @param hmm A list of parameters that specify the normal HMM, including
#'   `num_states`, `num_variables`, `num_subjects`, `mu`, `sigma`, `gamma`,
#'   `delta`.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
#'   treated as equal for all subjects.
#'
#' @return A list of matrices (one for each subject) of the forward variables.
#' @export
#'
#' @examples
norm_logforward <- function(x, hmm, state_dep_dist_pooled = FALSE) {
  n             <- nrow(x)
  num_states    <- hmm$num_states
  num_variables <- hmm$num_variables
  num_subjects  <- hmm$num_subjects
  lalpha        <- list()
  for (i in 1:num_subjects) {
    s_ind   <- i
    if (state_dep_dist_pooled) {
      s_ind <- 1
    }
    lalpha[[i]] <- matrix(NA, nrow = num_states, ncol = n)
    P           <- rep(1, num_states)
    for (j in 1:num_variables) {
      P <- P*stats::dnorm(x[1, j, i], hmm$mu[[j]][s_ind, ],
                          hmm$sigma[[j]][s_ind, ])
    }
    forward_probs     <- hmm$delta[[i]]*P
    sum_forward_probs <- sum(forward_probs)
    loglikelihood     <- log(sum_forward_probs)
    forward_probs     <- forward_probs/sum_forward_probs
    lalpha[[i]][, 1]  <- loglikelihood + log(forward_probs)

    for (t in 2:n) {
      P     <- rep(1, num_states)
      for (j in 1:num_variables) {
        if (!is.na(x[t, j, i])) {
          P <- P*stats::dnorm(x[t, j, i], hmm$mu[[j]][s_ind, ],
                              hmm$sigma[[j]][s_ind, ])
        }
      }
      forward_probs     <- forward_probs %*% hmm$gamma[[i]][[t]]*P
      sum_forward_probs <- sum(forward_probs)
      loglikelihood     <- loglikelihood + log(sum_forward_probs)
      forward_probs     <- forward_probs/sum_forward_probs
      lalpha[[i]][, t]  <- loglikelihood + log(forward_probs)
    }
  }
  lalpha
}
