#' Compute negative log-likelihood of normal HMM parameters
#'
#' This function computes the negative log-likelihood that the given normal
#' HMM parameters could have generated the data being fit.
#'
#' @param working_params A vector of the working normal parameters for the
#'   HMM.
#' @param x The data to be fit with an HMM in the form of a 3D array. The
#'   first index (row) corresponds to time, the second (column) to the
#'   variable number, and the third (matrix number) to the subject number.
#' @param design A list of design matrices for each subject with each row
#'   indicating the time and each column indicating the value of the
#'   covariate.
#' @param num_states The number of states in the desired HMM.
#' @param num_variables The number of variables in the data.
#' @param num_subjects The number of subjects that generated the data.
#' @param num_covariates The number of covariates in the data that the
#'   transition probability depends on.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
#'   treated as equal for all subjects.
#'
#' @return A number indicating the negative loglikelihood
#' @export

norm_loglikelihood <- function(working_params, x, design,
                               num_states, num_variables, num_subjects,
                               num_covariates,
                               state_dep_dist_pooled = FALSE) {
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
      P <- P*stats::dnorm(x[1, j, i], pn$mu[[j]][s_ind, ],
                          pn$sigma[[j]][s_ind, ])
    }
    forward_probs     <- pn$delta[[i]]*P
    sum_forward_probs <- sum(forward_probs)
    loglikelihood     <- log(sum_forward_probs)
    forward_probs     <- forward_probs/sum_forward_probs

    for (t in 2:num_time) {
      P     <- rep(1, num_states)
      for (j in 1:num_variables) {
        if (!is.na(x[t, j, i])) {
          P <- P*stats::dnorm(x[t, j, i], pn$mu[[j]][s_ind, ],
                              pn$sigma[[j]][s_ind, ])
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
