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
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
#'   treated as equal for all subjects.
#'
#' @return A list of the natural paramters.
#' @export
#' @examples
#' ### 1. Get the working parameters
#' # define values of parameters
#' num_states     <- 2
#' num_variables  <- 2
#' num_subjects   <- 2
#' num_covariates <- 2
#' mu             <- list(matrix(c(1, 5, 2, 4), 2, 2, byrow = TRUE),
#'                        matrix(c(1, 5, 2, 4), 2, 2, byrow = TRUE))
#' sigma          <- list(matrix(c(1, 2, 1, 1.5), 2, 2, byrow = TRUE),
#'                        matrix(c(1, 2, 1, 1.5), 2, 2, byrow = TRUE))
#' beta           <- matrix(c(-2, 0, 0), nrow = 1, ncol = 3)
#' delta          <- list(c(1/2, 1/2), c(1/2, 1/2))
#'
#' #transform to working parametes
#' working_params <- norm_working_params(num_states, num_variables,
#'                                       num_subjects, mu, sigma, beta, delta)
#' ### 2. Transform working_params to natural
#' norm_natural_params(num_states, num_variables, num_subjects, num_covariates,
#'                     working_params, state_dep_dist_pooled = FALSE)
norm_natural_params <- function(num_states, num_variables, num_subjects,
                                num_covariates, working_params,
                                state_dep_dist_pooled = FALSE) {
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
  }
  d <- split_vec(working_params, delta_start, delta_end, num_states - 1)
  for (i in 1:num_subjects) {
    foo        <- c(1, exp(d[[i]]))
    delta[[i]] <- foo/sum(foo)
  }
  list(mu = mu, sigma = sigma, beta = beta, delta = delta)
}
