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
#' @param state_dep_dist_pooled A logical variable indicating whether the
#'   state dependent distribution parameters `alpha` and `theta` should be
#'   treated as equal for all subjects.
#'
#' @return A list of the natural parameters.
#' @export
gam0_natural_params <- function(num_states, num_variables, num_subjects,
                                num_covariates, working_params,
                                state_dep_dist_pooled = FALSE) {
  ns             <- num_subjects
  alpha_start    <- 1
  alpha_end      <- num_states*num_variables*num_subjects
  theta_end      <- alpha_end + num_states*num_variables*num_subjects
  zweight_end    <- theta_end + num_subjects*num_variables
  alpha_len      <- theta_len <- num_subjects*num_states
  zweight_len    <- num_variables*num_subjects
  if (state_dep_dist_pooled) {
    alpha_end    <- num_states*num_variables
    theta_end    <- alpha_end + num_states*num_variables
    zweight_end  <- theta_end + num_variables
    alpha_len    <- theta_len <- num_states
    zweight_len  <- num_variables
  }
  theta_start   <- alpha_end + 1
  zweight_start <- theta_end + 1
  beta_start    <- zweight_end  + 1
  beta_end      <- zweight_end  +
    (num_states^2 - num_states)*(num_covariates + 1)
  delta_start   <- beta_end + 1
  delta_end     <- length(working_params)

  alpha   <- split_vec(working_params, alpha_start, alpha_end, alpha_len,
                       exp = TRUE)
  theta   <- split_vec(working_params, theta_start, theta_end, theta_len,
                       exp = TRUE)
  zweight <- split_vec(working_params, zweight_start, zweight_end, zweight_len,
                       tanh = TRUE)
  for (j in 1:num_variables) {
    alpha[[j]] <- matrix(alpha[[j]], ncol = num_states, byrow = TRUE)
    theta[[j]] <- matrix(theta[[j]], ncol = num_states, byrow = TRUE)
  }
  beta <- matrix(working_params[beta_start:beta_end],
                 nrow = num_states^2 - num_states)
  delta <- list()
  if (num_states == 1) {
    for (i in 1:num_subjects) {
      delta[[i]] = 1
      return(list(alpha = alpha, theta = theta, gamma = gamma, delta = delta))
    }
  }
  d <- split_vec(working_params, delta_start, delta_end, num_states - 1)
  for (i in 1:num_subjects) {
    foo        <- c(1, exp(d[[i]]))
    delta[[i]] <- foo/sum(foo)
  }
  list(alpha = alpha, theta = theta, zweight = zweight,
       beta = beta, delta = delta)
}
