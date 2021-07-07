#' Generate data from a normal HMM
#'
#' This function generates data from a normal HMM with specified parameters.
#'
#' @inheritParams norm_viterbi
#' @inheritParams norm_loglikelihood
#' @param num_sample The size of the desired sample (number of timesteps).
#'
#' @return A list of the data and the states that generated the data.
#' @export

norm_generate_sample <- function(num_sample, hmm, design,
                                 state_dep_dist_pooled = FALSE) {
  state_vec     <- 1:hmm$num_states
  num_variables <- hmm$num_variables
  num_subjects  <- hmm$num_subjects
  state         <- matrix(numeric(num_sample*num_subjects),
                          nrow = num_sample, ncol = num_subjects)
  gamma         <- norm_gamma(hmm$num_states, num_subjects, num_sample,
                              beta, design)
  for (i in 1:num_subjects) {
    state[1, i]   <- sample(state_vec, 1, prob = hmm$delta[[i]])
    for (t in 2:num_sample) {
      state[t, i] <- sample(state_vec, 1,
                            prob = hmm$gamma[[i]][[t]][state[(t - 1), i], ])
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
