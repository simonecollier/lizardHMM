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
#' @param pn A list containing the natural parameters for the gamma HMM.
#' @param state_dep_dist_pooled A logical variable indicating whether the
#'   state dependent distribution parameters `alpha` and `theta` should be
#'   treated as equal for all subjects.
#'
#' @return
#' @export
#'
#' @examples
gam0_allprobs <- function(num_states, num_variables, num_subjects, num_time,
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
          if (x[t, j, i] == 0) {
            P <- P*c(pn$zweight[[j]][s_ind], stats::dgamma(x[t, j, i],
                                            pn$alpha[[j]][s_ind, 2:num_states],
                                            pn$theta[[j]][s_ind, 2:num_states]))
          } else {
            P <- P*stats::dgamma(x[t, j, i], pn$alpha[[j]][s_ind, ],
                                 pn$theta[[j]][s_ind, ])*
              c(1 - pn$zweight[[j]][s_ind], 1, 1)
          }
        }
      }
      allprobs[[i]][t, ] <- P
    }
  }
  allprobs
}
