norm_loglikelihood <- function(working_params, x, design,
                               num_states, num_variables, num_subjects,
                               num_covariates,
                               state_dep_dist_pooled = FALSE) {
  #' Compute negative log-likelihood of normal HMM parameters
  #'
  #' This function computes the negative log-likelihood that the given normal
  #' HMM parameters could have generated the data being fit.
  #'
  #' @inheritParams norm_natural_params
  #' @param x The data to be fit with an HMM in the form of a 3D array. The
  #'   first index (row) corresponds to time, the second (column) to the
  #'   variable number, and the third (matrix number) to the subject number.
  #' @param design A list of design matrices for each subject with each row
  #'   indicating the time and each column indicating the value of the
  #'   covariate.
  #'
  #' @return A number indicating the negative loglikelihood
  #' @export

  num_time  <- nrow(x)
  pn        <- norm_natural_params(num_states, num_variables, num_subjects,
                                   num_covariates, working_params,
                                   state_dep_dist_pooled)
  gamma     <- norm_gamma(num_states, num_subjects, num_time, pn$beta, design)
  cum_loglikelihood <- 0
  for (i in 1:num_subjects) {
    s_ind   <- i
    if (state_dep_dist_pooled) {
      s_ind <- 1
    }
    P   <- rep(1, num_states)
    for (j in 1:num_variables) {
      P <- P*dnorm(x[1, j, i], pn$mu[[j]][s_ind, ], pn$sigma[[j]][s_ind, ])
    }
    forward_probs     <- pn$delta[[i]]*P
    sum_forward_probs <- sum(forward_probs)
    loglikelihood     <- log(sum_forward_probs)
    forward_probs     <- forward_probs/sum_forward_probs

    for (t in 2:num_time) {
      P     <- rep(1, num_states)
      for (j in 1:num_variables) {
        if (!is.na(x[t, j, i])) {
          P <- P*dnorm(x[t, j, i], pn$mu[[j]][s_ind, ], pn$sigma[[j]][s_ind, ])
        }
      }
      forward_probs     <- forward_probs %*% gamma[[i]][[t]]*P
      sum_forward_probs <- sum(forward_probs)
      loglikelihood     <- loglikelihood + log(sum_forward_probs)
      forward_probs     <- forward_probs/sum_forward_probs
    }
    cum_loglikelihood   <- cum_loglikelihood + loglikelihood
  }
  - cum_loglikelihood
}
