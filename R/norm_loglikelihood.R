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
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_covariates The number of covariates in the data that the
#'   transition probability matrix depends on.
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
  allprobs  <- norm_allprobs(num_states, num_variables, num_subjects, num_time,
                             x, pn, state_dep_dist_pooled)
  gamma     <- fit_tpm(num_states, num_subjects, num_covariates, num_time,
                       pn$beta, design)
  cum_loglikelihood <- 0
  for (i in num_subjects) {
    delta <- matrix(pn$delta[[i]], ncol = num_states)
    if (num_covariates == 0) {
      loglikelihood <- foralg(num_time, num_states,
                              delta, gamma[[i]], allprobs)
    } else {
      loglikelihood <- foralg_covar(num_time, num_states,
                                    delta, gamma[[i]], allprobs)
    }
    cum_loglikelihood <- cum_loglikelihood + loglikelihood
  }
  - cum_loglikelihood
}
