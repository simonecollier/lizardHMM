#' Generate data from a normal HMM
#'
#' This function generates data from a normal HMM with specified parameters.
#'
#' @param num_sample The size of the desired sample (number of timesteps).
#' @param hmm A list of parameters that specify the normal HMM, including
#'   `num_states`, `num_variables`, `num_subjects`, `num_covariates`, `mu`,
#'   `sigma`, `beta`, `delta`.
#' @param design A list of design matrices for each subject with each row
#'   indicating the time and each column indicating the value of the
#'   covariate.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
#'   treated as equal for all subjects.
#'
#' @return A list of the data and the states that generated the data.
#' @export
norm_generate_sample <- function(num_sample, hmm, design,
                                 state_dep_dist_pooled = FALSE) {
  state_vec      <- 1:hmm$num_states
  num_variables  <- hmm$num_variables
  num_subjects   <- hmm$num_subjects
  num_covariates <- hmm$num_covariates
  state          <- matrix(0, nrow = num_sample, ncol = num_subjects)
  gamma          <- fit_tpm(hmm$num_states, num_subjects, num_covariates, num_time,
                           beta, design)
  for (i in 1:num_subjects) {
    state[1, i]   <- sample(state_vec, 1, prob = hmm$delta[[i]])
    for (t in 2:num_sample) {
      if (num_covariates != 0) {
        state[t, i] <- sample(state_vec, 1,
                              prob = gamma[[i]][state[(t - 1), i], , t])
      } else {
        state[t, i] <- sample(state_vec, 1,
                              prob = gamma[[i]][state[(t - 1), i], ])
      }
    }
  }
  x <- array(dim = c(num_sample, num_variables, num_subjects))
  for (j in 1:num_variables) {
    for (i in 1:num_subjects) {
      s_ind   <- i
      if (state_dep_dist_pooled) {
        s_ind <- 1
      }
      x[, j, i] <- stats::rnorm(num_sample,
                                mean = hmm$mu[[j]][s_ind, state[, i]],
                                sd = hmm$sigma[[j]][s_ind, state[, i]])
    }
  }
  list(state = state, observ = x)
}
