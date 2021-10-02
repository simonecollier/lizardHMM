#' Fit an HMM
#'
#' This function fits data with an HMM by maximizing the likelihood estimate
#' given initial normal parameters.
#'
#' @param x The data to be fit with an HMM in the form of a 3D array. The
#'   first index (row) corresponds to time, the second (column) to the
#'   variable number, and the third (matrix number) to the subject number.
#' @param design A list of design matrices for each subject with each row
#'   indicating the time and each column indicating the value of the
#'   covariate.
#' @param num_states The number of states in the desired HMM.
#' @param num_variables The number of variables in the data.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_covariates The number of covariates in the data that the
#'   transition probability matrix depends on.
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
#' @param beta0 A matrix of the initial regression coefficients for the effect
#'    of the covariates on the transition probability matrices `gamma`.
#' @param delta0 A list with each element being the starting initial state
#'   distribution vector of the HMM for the subject corresponding to that
#'   index.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
#'   treated as equal for all subjects.
#' @param iterlim A value indicating the number of iterations `nlm()` should run
#'   before exiting.
#' @param hessian A logical variable indicating whether to compute the hessian
#'   at the minimum.
#'
#' @return A list of parameters that specify the fitted normal HMM.
#' @export
norm_fit_hmm <- function(x, design, num_states, num_variables, num_subjects,
                         num_covariates,
                         mu0, sigma0, beta0, delta0,
                         state_dep_dist_pooled = FALSE, iterlim = 100,
                         hessian = FALSE, print.level = 0) {
  num_time       <- nrow(x)
  working_params <- norm_working_params(num_states, num_variables, num_subjects,
                                        mu0, sigma0, beta0, delta0)
  hmm   <- stats::nlm(norm_loglikelihood,
                      working_params,
                      x = x,
                      design = design,
                      num_states = num_states,
                      num_variables = num_variables,
                      num_subjects = num_subjects,
                      num_covariates = num_covariates,
                      state_dep_dist_pooled = state_dep_dist_pooled,
                      iterlim = iterlim,
                      hessian = hessian,
                      print.level = print.level)
  pn    <- norm_natural_params(num_states = num_states,
                               num_variables = num_variables,
                               num_subjects = num_subjects,
                               num_covariates = num_covariates,
                               working_params = hmm$estimate,
                               state_dep_dist_pooled = state_dep_dist_pooled)
  gamma <- fit_tpm(num_states, num_subjects, num_covariates, num_time,
                   pn$beta, design)

  mllk <- hmm$minimum
  p    <- length(working_params)
  n    <- sum(!is.na(x))
  AIC  <- 2*(mllk + p)
  BIC  <- 2*mllk + p*log(n)

  if (hessian) {
    h <- hmm$hessian
    if (num_states != 1) {
      if (state_dep_dist_pooled) {
        d <- num_states - 1
      } else {
        d <- (num_states - 1)*num_subjects
      }
      h <- h[1:(nrow(h) - d), 1:(nrow(h) - d)]
    }
    if (det(h) != 0) {
      h <- solve(h)
      return(list(num_states = num_states,
                  num_variables = num_variables,
                  num_subjects = num_subjects,
                  num_covariates = num_covariates,
                  mu = pn$mu,
                  sigma = pn$sigma,
                  beta = pn$beta,
                  delta = pn$delta,
                  gamma = gamma,
                  working_params = hmm$estimate,
                  code = hmm$code,
                  max_loglikelihood = hmm$minimum,
                  AIC = AIC,
                  BIC = BIC,
                  inverse_hessian = h))
    }
    return(list(num_states = num_states,
                num_variables = num_variables,
                num_subjects = num_subjects,
                num_covariates = num_covariates,
                mu = pn$mu,
                sigma = pn$sigma,
                beta = pn$beta,
                delta = pn$delta,
                gamma = gamma,
                working_params = hmm$estimate,
                code = hmm$code,
                max_loglikelihood = hmm$minimum,
                AIC = AIC,
                BIC = BIC,
                hessian = hmm$hessian))
  }
  list(num_states = num_states,
       num_variables = num_variables,
       num_subjects = num_subjects,
       num_covariates = num_covariates,
       mu = pn$mu,
       sigma = pn$sigma,
       delta = pn$delta,
       beta = pn$beta,
       gamma = gamma,
       working_params = hmm$estimate,
       code = hmm$code,
       max_loglikelihood = hmm$minimum,
       AIC = AIC,
       BIC = BIC)
}
