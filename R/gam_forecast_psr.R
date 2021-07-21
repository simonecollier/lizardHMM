#' Compute forecast gamma pseuso-residuals
#'
#' This function computes the gamma forecast pseudo-residuals of the data `x`
#'   fitted with `hmm` (with gamma state dependent distributions).
#'
#' @param x The data to be fit with an HMM in the form of a 3D array. The
#'   first index (row) corresponds to time, the second (column) to the
#'   variable number, and the third (matrix number) to the subject number.
#' @param hmm A list of parameters that specify the gamma HMM, including
#'   `num_states`, `num_variables`, `num_subjects`, `alpha`, `theta`, `gamma`,
#'   `delta`.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `alpha` and `theta` should be
#'   treated as equal for all subjects.
#'
#' @return A list of vectors (one for each subject) of the pseudo-residuals.
#' @export

gam_forecast_psr <- function(x, hmm, state_dep_dist_pooled = FALSE) {
  n             <- nrow(x)
  num_states    <- hmm$num_states
  num_variables <- hmm$num_variables
  num_subjects  <- hmm$num_subjects
  la            <- gam_logforward(x, hmm, state_dep_dist_pooled)
  forecast_psr  <- list()
  for (i in 1:num_subjects) {
    s_ind   <- i
    if (state_dep_dist_pooled) {
      s_ind <- 1
    }
    pstepmat          <- matrix(NA, n, num_states)
    forecast_psr[[i]] <- rep(NA, n)
    ind_step          <- 1:n
    for (j in 1:num_variables) {
      ind_step <- sort(intersect(ind_step, which(!is.na(x[, j, i]))))
    }
    for (k in ind_step) {
      for (m in 1:num_states) {
        P   <- 1
        for (j in 1:num_variables) {
          P <- P*stats::pgamma(x[k, j, i], shape = hmm$alpha[[j]][s_ind, m],
                               scale = hmm$theta[[j]][s_ind, m])
        }
        pstepmat[k, m] <- P
      }
    }
    if (1 %in% ind_step) {
      forecast_psr[[i]][1] <- stats::qnorm(hmm$delta[[i]] %*% pstepmat[1, ])
    }
    for (t in 2:n) {
      c <- max(la[[i]][, t - 1])
      a <- exp(la[[i]][, t - 1] - c)
      if (t %in% ind_step) {
        forecast_psr[[i]][t] <- stats::qnorm(t(a) %*%
                                               (hmm$gamma[[i]][[t]]/sum(a)) %*%
                                               pstepmat[t, ])
      }
    }
  }
  forecast_psr
}
