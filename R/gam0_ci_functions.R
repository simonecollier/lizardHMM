#' Convert working parameters to a vector of natural parameters
#'
#' This funciton converts the vecgtor of working parameters to a vector of
#' natural parameters not inculding delta.
#'
#' @param num_states The number of states in the desired HMM.
#' @param num_variables The number of variables in the data.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_covariates The number of covariates in the data that the
#'   transition probability matrix depends on.
#' @param working_params A vector of the working gamma parameters for the
#'   HMM.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `alpha` and `theta` should be
#'   treated as equal for all subjects.
#'
#' @return A vector of the natural parameters.
#' @export
gam0_natural_vec <- function(num_states, num_variables, num_subjects,
                            num_covariates, working_params,
                            state_dep_dist_pooled = FALSE) {
  alpha_start <- 1
  alpha_end   <- num_states*num_variables*num_subjects
  theta_end   <- alpha_end + num_states*num_variables*num_subjects
  zweight_end <- theta_end + num_subjects*num_variables
  alpha_len   <- theta_len <- num_subjects*num_states
  zweight_len <- num_variables*num_subjects
  if (state_dep_dist_pooled) {
    alpha_end   <- num_states*num_variables
    theta_end   <- alpha_end + num_states*num_variables
    zweight_end <- theta_end + num_variables
    alpha_len   <- theta_len <- num_states
    zweight_len  <- num_variables
  }
  theta_start   <- alpha_end + 1
  zweight_start <- theta_end + 1
  beta_start    <- zweight_end + 1
  beta_end     <- zweight_end + (num_states^2 - num_states)*(num_covariates + 1)

  natural <- numeric()
  natural[alpha_start:alpha_end] <- exp(working_params[alpha_start:alpha_end])
  natural[theta_start:theta_end] <- exp(working_params[theta_start:theta_end])
  natural[zweight_start:zweight_end] <- tanh(working_params[zweight_start:
                                                              zweight_end])
  natural[beta_start:beta_end] <- working_params[beta_start:beta_end]

  natural
}


#' Reformat confidence interval data
#'
#' This is a helper function for `gam_ci()` which reformats the output so that
#' it is more easily interpreted.
#'
#' @param num_states The number of states in the desired HMM.
#' @param num_variables The number of variables in the data.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_covariates The number of covariates in the data that the
#'   transition probability depends on.
#' @param estimate_vec A vector containing the estimated natural parameters
#'   of the gamma HMM in the format outputte by `gam_natural_vec()`.
#' @param upper_vec A vector containing the upper confidence interval of the
#'   estimated natural parameters of the gamma HMM in the format outputte by
#'   `gam_natural_vec()`.
#' @param lower_vec A vector containing the lower confidence interval of the
#'   estimated natural parameters of the gamma HMM in the format outputte by
#'   `gam_natural_vec()`.
#'@param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `alpha` and `theta` should be
#'   treated as equal for all subjects.
#'
#' @return A list containing the upper and lower confidence interval for each
#'   estimated parameter.
#' @export
gam0_ci_data <- function(num_states, num_variables, num_subjects,
                         num_covariates, estimate_vec, upper_vec, lower_vec,
                         state_dep_dist_pooled = FALSE) {
  alpha_start <- 1
  alpha_end   <- num_states*num_variables*num_subjects
  theta_end   <- alpha_end + num_states*num_variables*num_subjects
  zweight_end <- theta_end + num_subjects*num_variables
  alpha_len   <- theta_len <- num_subjects*num_states
  zweight_len <- num_variables*num_subjects
  if (state_dep_dist_pooled) {
    alpha_end   <- num_states*num_variables
    theta_end   <- alpha_end + num_states*num_variables
    zweight_end <- theta_end + num_variables
    alpha_len   <- theta_len <- num_states
    zweight_len  <- num_variables
  }
  theta_start   <- alpha_end + 1
  zweight_start <- theta_end + 1
  beta_start    <- zweight_end + 1
  beta_end      <- zweight_end + (num_states^2 - num_states)*(num_covariates + 1)

  alpha_estimate <- split_vec(estimate_vec, alpha_start, alpha_end, alpha_len)
  alpha_upper    <- split_vec(upper_vec, alpha_start, alpha_end, alpha_len)
  alpha_lower    <- split_vec(lower_vec, alpha_start, alpha_end, alpha_len)
  theta_estimate <- split_vec(estimate_vec, theta_start, theta_end, theta_len)
  theta_upper    <- split_vec(upper_vec, theta_start, theta_end, theta_len)
  theta_lower    <- split_vec(lower_vec, theta_start, theta_end, theta_len)
  zweight_estimate <- split_vec(estimate_vec, zweight_start, zweight_end,
                                zweight_len)
  zweight_upper    <- split_vec(upper_vec, zweight_start, zweight_end,
                                zweight_len)
  zweight_lower    <- split_vec(zweight_vec, zweight_start, zweight_end,
                                zweight_len)
  for (j in 1:num_variables) {
    alpha_estimate[[j]] <- matrix(alpha_estimate[[j]], ncol = num_states,
                                  byrow = TRUE)
    alpha_upper[[j]]    <- matrix(alpha_upper[[j]], ncol = num_states,
                                  byrow = TRUE)
    alpha_lower[[j]]    <- matrix(alpha_lower[[j]], ncol = num_states,
                                  byrow = TRUE)
    theta_estimate[[j]] <- matrix(theta_estimate[[j]], ncol = num_states,
                                  byrow = TRUE)
    theta_upper[[j]]    <- matrix(theta_upper[[j]], ncol = num_states,
                                  byrow = TRUE)
    theta_lower[[j]]    <- matrix(theta_lower[[j]], ncol = num_states,
                                  byrow = TRUE)
  }
  alpha         <- list(estimate = alpha_estimate,
                        upper = alpha_upper,
                        lower = alpha_lower)
  theta         <- list(estimate = theta_estimate,
                        upper = theta_upper,
                        lower = theta_lower)
  zweight       <- list(estimate = zweight_estimate,
                        upper = zweight_upper,
                        lower = zweight_lower)

  beta_estimate <- matrix(estimate_vec[beta_start:beta_end],
                          nrow = num_states^2 - num_states)
  beta_upper    <- matrix(upper_vec[beta_start:beta_end],
                          nrow = num_states^2 - num_states)
  beta_lower    <- matrix(lower_vec[beta_start:beta_end],
                          nrow = num_states^2 - num_states)
  beta          <- list(estimate = beta_estimate,
                        upper = beta_upper,
                        lower = beta_lower)

  list(alpha = alpha, theta = theta, zweight = zweight, beta = beta)
}
