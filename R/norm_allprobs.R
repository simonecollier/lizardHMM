#' Compute allprobs
#'
#' This function computes the probability of each observation arising from each
#' state.
#'
#' @param num_states The number of states in the desired HMM.
#' @param num_variables The number of variables in the data.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_time A value indicating the length of the data.
#' @param x The data to be fit with an HMM in the form of a 3D array. The
#'   first index (row) corresponds to time, the second (column) to the
#'   variable number, and the third (matrix number) to the subject number.
#' @param pn A list conatining the natural parameters for the normal HMM.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
#'   treated as equal for all subjects.
#'
#' @return
#' @export
#'
#' @examples
norm_allprobs <- function(num_states, num_variables, num_subjects, num_time,
                          x, pn, state_dep_dist_pooled = FALSE) {
  allprobs   <- list()
  for (i in 1:num_subjects) {
    s_ind    <- i
    if (state_dep_dist_pooled) {
      s_ind  <- 1
    }
    allprobs[[i]] <- matrix(0, ncol = num_states, nrow = num_time)
    for (t in 1:num_time) {
      P     <- rep(1, num_states)
      for (j in 1:num_variables) {
        if (!is.na(x[t, j, i])) {
          P <- P*stats::dnorm(x[t, j, i], pn$mu[[j]][s_ind, ],
                              pn$sigma[[j]][s_ind, ])
        }
      }
      allprobs[[i]][t, ] <- P
    }
  }
  allprobs
}
