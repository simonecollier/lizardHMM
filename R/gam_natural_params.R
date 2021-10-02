#' Transform gamma HMM parameters from working to natural
#'
#' This function transforms the working gamma HMM parameters back into the
#' original format of the natural parameters and outputs them as a list. This
#' function is the reverse of `norm_working_params()`.
#'
#' @param num_states The number of states in the desired HMM.
#' @param num_variables The number of variables in the data.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_covariates The number of covariates in the data that the
#'   transition probability matrix depends on.
#' @param working_params A vector of the working gamma parameters for the
#'   HMM as outputted by `norm_working_params()`.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `alpha` and `theta` should be
#'   treated as equal for all subjects.
#'
#' @return A list of the natural paramters.
#' @export
gam_natural_params <- function(num_states, num_variables, num_subjects,
                               num_covariates, working_params,
                               state_dep_dist_pooled = FALSE) {

  ind <- gam_working_ind(num_states, num_variables, num_subjects,
                         num_covariates, state_dep_dist_pooled = FALSE)

  alpha <- split_vec(working_params, ind$alpha_start, ind$alpha_end,
                     ind$alpha_len, exp = TRUE)
  theta <- split_vec(working_params, ind$theta_start, ind$theta_end,
                     ind$theta_len, exp = TRUE)
  for (j in 1:num_variables) {
    alpha[[j]] <- matrix(alpha[[j]], ncol = num_states, byrow = TRUE)
    theta[[j]] <- matrix(theta[[j]], ncol = num_states, byrow = TRUE)
  }
  beta <- matrix(working_params[ind$beta_start:ind$beta_end],
                 nrow = num_states^2 - num_states)
  delta <- list()
  if (num_states == 1) {
    for (i in 1:num_subjects) {
      delta[[i]] = 1
      return(list(alpha = alpha, theta = theta, gamma = gamma, delta = delta))
    }
  }
  d <- split_vec(working_params, ind$delta_start, ind$delta_end, num_states - 1)
  for (i in 1:num_subjects) {
    foo        <- c(1, exp(d[[i]]))
    delta[[i]] <- foo/sum(foo)
  }
  list(alpha = alpha, theta = theta, beta = beta, delta = delta)
}
