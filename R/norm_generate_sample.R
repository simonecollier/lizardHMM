#' Generate data from a normal HMM
#'
#' This function generates data from a normal HMM with specified parameters.
#'
#' @param num_sample The size of the desired sample (number of timesteps).
#' @param hmm A list of parameters that specify the normal HMM, including
#'   `num_states`, `num_variables`, `num_subjects`, `mu`, `sigma`, `beta`,
#'   `delta`.
#' @param design A list of design matrices for each subject with each row
#'   indicating the time and each column indicating the value of the
#'   covariate.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
#'   treated as equal for all subjects.
#'
#' @return A list of the data and the states that generated the data.
#' @export
#' @examples
#' # define parameters
#' hmm <-  list(num_states = 3,
#'              num_variables = 2,
#'              num_subjects = 1,
#'              mu = list(matrix(c(12, 18, 22),
#'                        ncol = 3, nrow = 1, byrow = TRUE),
#'                        matrix(c(-12, -7, 0),
#'                               ncol = 3, nrow = 1, byrow = TRUE)),
#'              sigma = list(matrix(c(3, 1, 1.5),
#'                                  ncol = 3, nrow = 1, byrow = TRUE),
#'                           matrix(c(1, 3, 2),
#'                                  ncol = 3, nrow = 1, byrow = TRUE)),
#'              beta  = matrix(c(0.01, 0.02, 0.001,
#'                               0.01, 0.03, 0.004,
#'                               0.01, 0.01, 0.003,
#'                               0.01, 0.04, 0.002,
#'                               0.01, 0.01, 0.004,
#'                               0.01, 0.03, 0.001),
#'                               ncol = 3, nrow = 6, byrow = TRUE),
#'              delta = list(c(0.3, 0.2, 0.5)))
#'
#' num_sample       <- 1000
#' design           <- list(matrix(0, nrow = 1000, ncol = 3))
#' design[[1]][, 1] <- 1 # First column is the intercept
#' design[[1]][, 2] <- sample(c(0, 1), size = 1000,
#'                            prob = c(0.3, 0.7), replace = TRUE)
#' design[[1]][, 3] <- rnorm(1000, mean = 5, sd = 1)
#'
#' # generate sample
#' norm_generate_sample(num_sample, hmm, design, state_dep_dist_pooled = FALSE)
norm_generate_sample <- function(num_sample, hmm, design,
                                 state_dep_dist_pooled = FALSE) {
  state_vec     <- 1:hmm$num_states
  num_variables <- hmm$num_variables
  num_subjects  <- hmm$num_subjects
  state         <- matrix(numeric(num_sample*num_subjects),
                          nrow = num_sample, ncol = num_subjects)
  gamma         <- fit_tpm(hmm$num_states, num_subjects, num_sample,
                              hmm$beta, design)
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
