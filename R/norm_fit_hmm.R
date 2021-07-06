norm_fit_hmm <- function(x, design, num_states, num_variables, num_subjects,
                         num_covariates,
                         mu0, sigma0, beta0, delta0,
                         state_dep_dist_pooled = FALSE) {
  #' Fit an HMM
  #'
  #' This function fits data with an HMM by maximizing the likelihood estimate
  #' given initial normal parameters.
  #'
  #' @inheritParams norm_loglikelihood
  #' @param mu0 The starting values for the means of the normally distributed
  #'   state dependent distributions of the HMM. `mu0` is a list of matrices,
  #'   each matrix corresponding to a different variable in the data being fit.
  #'   The columns of the matrices correspond to the state number and the rows
  #'   correspond to the subject number.
  #' @param sigma0 The starting values for the standard deviations of the
  #'   normally distributed state dependent distributions of the HMM. `sigma0`
  #'   is a list of matrices, each matrix corresponding to a different variable
  #'   in the data being fit. The columns of the matrices correspond to the
  #'   state number and the rows correspond to the subject number.
  #' @param beta0 A matrix of regression coefficients.
  #' @param delta0 A list with each element being the starting initial state
  #'   distribution vector of the HMM for the subject corresponding to that
  #'   index.
  #'
  #' @return A list of parameters that specify the fitted normal HMM, including
  #'   `num_states`, `num_variables`, `num_subjects`, `mu`, `sigma`, `gamma`,
  #'   `delta`.
  #' @export

  num_time <- nrow(x)
  working_params <- norm_working_params(num_states, num_subjects,
                                        mu0, sigma0, beta0, delta0)
  hmm <- stats::nlm(norm_loglikelihood,
                    working_params,
                    x = x,
                    design = design,
                    num_states = num_states,
                    num_variables = num_variables,
                    num_subjects = num_subjects,
                    num_covariates = num_covariates,
                    state_dep_dist_pooled = state_dep_dist_pooled)
                    # hessian = hessian
  pn    <- norm_natural_params(num_states = num_states,
                               num_variables = num_variables,
                               num_subjects = num_subjects,
                               num_covariates = num_covariates,
                               working_params = hmm$estimate,
                               state_dep_dist_pooled = state_dep_dist_pooled)
  gamma <- norm_gamma(num_states, num_subjects, num_time, pn$beta, design)

  # if (hessian) {
  #   h <- hmm$hessian
  #   if (det(h) != 0) {
  #     h <- solve(h)
  #     jacobian <- norm_jacobian(num_states, num_variables, num_subjects,
  #                               hmm$estimate, pn,
  #                               state_dep_dist_pooled = state_dep_dist_pooled)
  #     h <- t(jacobian)%*%h%*%jacobian
  #   }
  # }

  # if ((hessian) & (det(hmm$hessian) != 0)) {
  #   list(num_states = num_states,
  #        num_variables = num_variables,
  #        num_subjects = num_subjects,
  #        mu = pn$mu,
  #        sigma = pn$sigma,
  #        gamma = gamma,
  #        delta = pn$delta,
  #        code = hmm$code,
  #        max_loglikelihood = hmm$minimum,
  #        inverse_hessian = h)
  # } else if (hessian) {
  #   list(num_states = num_states,
  #        num_variables = num_variables,
  #        num_subjects = num_subjects,
  #        mu = pn$mu,
  #        sigma = pn$sigma,
  #        gamma = gamma,
  #        delta = pn$delta,
  #        code = hmm$code,
  #        max_loglikelihood = hmm$minimum,
  #        hessian = hmm$hessian,
  #        determinant = det(hmm$hessian))
  # }
  list(num_states = num_states,
       num_variables = num_variables,
       num_subjects = num_subjects,
       mu = pn$mu,
       sigma = pn$sigma,
       gamma = gamma,
       delta = pn$delta,
       beta = pn$beta,
       code = hmm$code,
       max_loglikelihood = hmm$minimum)
}
