#' Compute log forward probabilities
#'
#' This function computes the log forward probabilities of the data based on
#' the HMM hmm (with normal state dependent distributions).
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
norm_logforward <- function(x, hmm, state_dep_dist_pooled = FALSE) {
  num_time      <- nrow(x)
  num_states    <- hmm$num_states
  num_variables <- hmm$num_variables
  num_subjects  <- hmm$num_subjects
  la            <- list()
  allprobs      <- norm_allprobs(num_states, num_variables,
                                 num_subjects, num_time,
                                 x, hmm, state_dep_dist_pooled = FALSE)
  for (i in 1:num_subjects) {
    la[[i]] <- matrix(NA, nrow = num_states, ncol = num_time)

    forward_probs     <- hmm$delta[[i]]*allprobs[[i]][1, ]
    sum_forward_probs <- sum(forward_probs)
    loglikelihood     <- log(sum_forward_probs)
    forward_probs     <- forward_probs/sum_forward_probs
    la[[i]][, 1]      <- loglikelihood + log(forward_probs)

    for (t in 2:num_time) {
      if (num_covariates != 0) {
        forward_probs <- forward_probs %*%
                           hmm$gamma[[i]][, , t]*allprobs[[i]][t, ]
      } else {
        forward_probs <- forward_probs %*% hmm$gamma[[i]]*allprobs[[i]][t, ]
      }
      sum_forward_probs <- sum(forward_probs)
      loglikelihood     <- loglikelihood + log(sum_forward_probs)
      forward_probs     <- forward_probs/sum_forward_probs
      la[[i]][, t]      <- loglikelihood + log(forward_probs)
    }
  }
  la
}
