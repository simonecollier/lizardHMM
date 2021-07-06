norm_natural_params <- function(num_states, num_variables, num_subjects,
                                num_covariates, working_params,
                                state_dep_dist_pooled = FALSE) {
  #' Transform normal HMM parameters from working to natural
  #'
  #' This function transforms the working normal HMM parameters back into the
  #' original format of the natural parameters and outputs them as a list. This
  #' function is the reverse of `norm_working_params()`.
  #'
  #' @inheritParams norm_working_params
  #' @param num_covariates The number of covariates in the data that the
  #'   transition probability depends on.
  #' @param working_params A vector of the working normal parameters for the
  #'   HMM.
  #' @param state_dep_dist_pooled A logical variable indiacting whether the
  #'   state dependent distribution parameters `mu` and `sigma` should be
  #'   treated as equal for all subjects.
  #'
  #' @return A list of the natural paramters.
  #' @export
  #' @examples
  #' norm_natural_params(2, 2, 2, 2, c(1, 5, 2, 4, 1, 5, 2, 4, 0, 0.6931472, 0,
  #'   0.4054651, 0, 0.6931472, 0, 0.4054651, -2, 0, 0, 0, 0))

  ns          <- num_subjects
  mu_start    <- 1
  mu_end      <- num_states*num_variables*num_subjects
  sigma_end   <- mu_end + num_states*num_variables*num_subjects
  mu_len      <- sigma_len <- num_subjects*num_states
  if (state_dep_dist_pooled) {
    mu_end    <- num_states*num_variables
    sigma_end <- mu_end + num_states*num_variables
    mu_len    <- sigma_len <- num_states
  }
  sigma_start <- mu_end + 1
  beta_start  <- sigma_end + 1
  beta_end    <- sigma_end + (num_states^2 - num_states)*(num_covariates + 1)
  delta_start <- beta_end + 1
  delta_end   <- length(working_params)

  mu    <- split_vec(working_params, mu_start, mu_end, mu_len)
  sigma <- split_vec(working_params, sigma_start, sigma_end, sigma_len,
                     exp = TRUE)
  for (j in 1:num_variables) {
    mu[[j]]    <- matrix(mu[[j]], ncol = num_states, byrow = TRUE)
    sigma[[j]] <- matrix(sigma[[j]], ncol = num_states, byrow = TRUE)
  }
  beta <- matrix(working_params[beta_start:beta_end],
                 nrow = num_states^2 - num_states)
  delta <- list()
  if (num_states == 1) {
    for (i in 1:num_subjects) {
      delta[[i]] = 1
      return(list(mu = mu, sigma = sigma, gamma = gamma, delta = delta))
    }
    d <- split_vec(working_params, delta_start, delta_end, num_states - 1)
    for (i in 1:num_subjects) {
      foo        <- c(1, exp(d[[i]]))
      delta[[i]] <- foo/sum(foo)
    }
  }
  list(mu = mu, sigma = sigma, beta = beta, delta = delta)
}
