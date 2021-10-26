#' Get indices of parameters in working vector
#'
#' This function finds the indices for each parameter of the normal HMM within
#' the vector of working parameters outputted by `norm_working_params()`.
#'
#' @param num_states The number of states in the desired HMM.
#' @param num_variables The number of variables in the data.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_covariates The number of covariates in the data that the
#'   transition probability matrix depends on.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
#'   treated as equal for all subjects.
#'
#' @return A list of the start and end indices.
#' @export
norm0_working_ind <- function(num_states, num_variables, num_subjects,
                              num_covariates, state_dep_dist_pooled = FALSE) {
  mu_start    <- 1
  mu_end      <- num_states*num_variables*num_subjects
  sigma_end   <- mu_end + num_states*num_variables*num_subjects
  zweight_end <- sigma_end + num_subjects*num_variables
  mu_len      <- sigma_len <- num_subjects*num_states
  zweight_len <- num_variables*num_subjects
  if (state_dep_dist_pooled) {
    mu_end      <- num_states*num_variables
    sigma_end   <- mu_end + num_states*num_variables
    zweight_end <- sigma_end + num_variables
    mu_len      <- sigma_len <- num_states
    zweight_len <- num_variables
  }
  sigma_start   <- mu_end + 1
  zweight_start <- sigma_end + 1
  beta_start    <- zweight_end + 1
  beta_end     <- zweight_end + (num_states^2 - num_states)*(num_covariates + 1)
  delta_start  <- beta_end + 1
  delta_end    <- delta_start + num_states - 2

  list(mu_start    = mu_start,
       mu_end      = mu_end,
       mu_len      = mu_len,
       sigma_start = sigma_start,
       sigma_end   = sigma_end,
       sigma_len   = sigma_len,
       beta_start  = beta_start,
       beta_end    = beta_end,
       delta_start = delta_start,
       delta_end   = delta_end)
}


#' Transform normal HMM parameters from working to natural
#'
#' This function transforms the working normal HMM parameters back into the
#' original format of the natural parameters and outputs them as a list. This
#' function is the reverse of `norm_working_params()`.
#'
#' @param num_states The number of states in the desired HMM.
#' @param num_variables The number of variables in the data.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_covariates The number of covariates in the data that the
#'   transition probability matrix depends on.
#' @param working_params A vector of the working normal parameters for the
#'   HMM as outputted by `norm_working_params()`.
#' @param state_dep_dist_pooled A logical variable indicating whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
#'   treated as equal for all subjects.
#'
#' @return A list of the natural parameters.
#' @export
norm0_natural_params <- function(num_states, num_variables, num_subjects,
                                 num_covariates, working_params,
                                 state_dep_dist_pooled = FALSE) {

  ind <- norm_working_ind(num_states, num_variables, num_subjects,
                          num_covariates, state_dep_dist_pooled = FALSE)

  mu    <- split_vec(working_params, ind$mu_start, ind$mu_end, ind$mu_len)
  sigma <- split_vec(working_params, ind$sigma_start, ind$sigma_end,
                     ind$sigma_len, exp = TRUE)
  for (j in 1:num_variables) {
    mu[[j]]    <- matrix(mu[[j]], ncol = num_states, byrow = TRUE)
    sigma[[j]] <- matrix(sigma[[j]], ncol = num_states, byrow = TRUE)
  }
  beta <- matrix(working_params[ind$beta_start:ind$beta_end],
                 nrow = num_states^2 - num_states)
  delta <- list()
  if (num_states == 1) {
    for (i in 1:num_subjects) {
      delta[[i]] = 1
      return(list(mu = mu, sigma = sigma, gamma = gamma, delta = delta))
    }
  }
  d <- split_vec(working_params, ind$delta_start, ind$delta_end, num_states - 1)
  for (i in 1:num_subjects) {
    foo        <- c(1, exp(d[[i]]))
    delta[[i]] <- foo/sum(foo)
  }
  list(mu = mu, sigma = sigma, beta = beta, delta = delta)
}
