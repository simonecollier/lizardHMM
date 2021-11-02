#' Compute transition probability matrices for plotting
#'
#' This function computes the transition probability matrix entries with their
#' confidence intervals for different covariate values in order to plot how
#' they change. THIS FUNCTION DOES NOT SUPPORT UNPOOLED TRANSITION PROBABILITY
#' MATRICES.
#'
#' @param num_states The number of states in the desired HMM.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_covariates The number of covariates in the data that the
#'   transition probability matrix depends on.
#' @param len_covariates A value indicating how many covariate values will be
#'   used.
#' @param design A list of design matrices, one for each subject which
#'   indicate the values of the each of the covariates (column) at each point
#'   in time (row).
#' @param conf_intervals A list of the confidence intervals for each fitted
#'   parameter of the HMM as outputted by the functions `norm_ci`, `gam_ci`, or
#'   `gam0_ci`.
#'
#' @return A list of 3 matrices, one for the estimate and the upper, and lower
#'   confidence interval, with the rows of the matrices containing the entries
#'   of the transition probability matrices at each covariate value.
#' @export
covariate_ci <- function(hmm, len_covariates, design, n = 100, level = 0.975,
                         state_dep_dist_pooled = FALSE) {

  inds <- gam0_working_ind(hmm$num_states, hmm$num_variables, hmm$num_subjects,
                           hmm$num_covariates, state_dep_dist_pooled)
  beta_sd       <- sqrt(diag(hmm$inverse_hessian))[inds$beta_start:inds$beta_end]
  beta_vec      <- hmm$working_params[inds$beta_start:inds$beta_end]
  len           <- length(beta_sd)
  beta_entries  <- matrix(0, ncol = len, nrow = n)
  for (l in 1:n) {
    sample            <- rnorm(len, mean = beta_vec, sd = beta_sd)
    beta_entries[l, ] <- sample
  }
  gamma_entries <- list()
  delta_entries <- list()
  for (t in 1:len_covariates) {
    gamma_entries[[t]] <- matrix(0, ncol = hmm$num_states*hmm$num_states,
                                 nrow = n)
    delta_entries[[t]] <- matrix(0, ncol = hmm$num_states, nrow = n)
    for (l in 1:n) {
      beta <- matrix(beta_entries[l, ], nrow = num_states^2 - num_states)
      gamma <- fit_tpm(hmm$num_states, hmm$num_subjects, hmm$num_covariates,
                       1, beta, list(matrix(design[[1]][t, ], nrow = 1)))
      gamma_entries[[t]][l, ] <- unlist(gamma)
      eigs <- eigen(x = t(gamma[[1]][, , 1]))
      values <- numeric()
      for (value in eigs$value) {
        if (zapsmall(Im(value)) == 0) {
          values <- c(values, zapsmall(Re(value)))
        } else {
          values <- c(values, zapsmall(value))
        }
      }
      ind <- which(values == 1)
      delta_entries[[t]][l, ] <- Re(eigs$vectors[, ind]/sum(eigs$vectors[, ind]))
    }
  }
  upper_gamma <- matrix(0, ncol = hmm$num_states*hmm$num_states,
                        nrow = len_covariates)
  lower_gamma <- matrix(0, ncol = hmm$num_states*hmm$num_states,
                        nrow = len_covariates)
  for (i in 1:(hmm$num_states*hmm$num_states)) {
    for (t in 1:len_covariates) {
      upper_gamma[t, i] <- quantile(gamma_entries[[t]][, i],
                                    probs = level, na.rm = TRUE)
      lower_gamma[t, i] <- quantile(gamma_entries[[t]][, i],
                                    probs = 1 - level, na.rm = TRUE)
    }
  }
  upper_delta <- matrix(0, ncol = hmm$num_states, nrow = len_covariates)
  lower_delta <- matrix(0, ncol = hmm$num_states, nrow = len_covariates)
  for (j in 1:num_states) {
    for (t in 1:len_covariates) {
      upper_delta[t, j] <- quantile(delta_entries[[t]][, j],
                                    probs = level, na.rm = TRUE)
      lower_delta[t, j] <- quantile(delta_entries[[t]][, j],
                                    probs = 1 - level, na.rm = TRUE)
    }
  }
  gamma <- fit_tpm(hmm$num_states, hmm$num_subjects, hmm$num_covariates,
                   len_covariates, hmm$beta, design)
  gamma_mat <- matrix(unlist(gamma), ncol = hmm$num_states*hmm$num_states,
                             byrow = TRUE)
  delta_mat <- matrix(0, ncol = hmm$num_states, nrow = len_covariates)
  for (t in 1:len_covariates){
    eigs <- eigen(t(gamma[[1]][, , t]))
    values <- numeric()
    for (value in eigs$value) {
      if (zapsmall(Im(value)) == 0) {
        values <- c(values, zapsmall(Re(value)))
      } else {
        values <- c(values, zapsmall(value))
      }
    }
    ind <- which(values == 1)
    delta_mat[t, ] <- Re(eigs$vectors[, ind]/sum(eigs$vectors[, ind]))
  }

  list(gammas = gamma_mat,
       upper_gamma = upper_gamma,
       lower_gamma = lower_gamma,
       deltas = delta_mat,
       upper_delta = upper_delta,
       lower_delta = lower_delta)
}


#' Plot entries of transition probability matrices
#'
#' This function plots the transition probability matrix entries with their
#' confidence intervals for different covariate values.
#'
#' @param num_states The number of states in the HMM.
#' @param covar_vec The vector of covariate values to plot the entries over.
#' @param tpm_entries_list A list of 3 matrices, one for the estimate and the
#'   upper, and lower confidence interval, with the rows of the matrices
#'   containing the entries of the transition probability matrices at each
#'   covariate value.
#' @param covariate_name Label for the x-axis.
#'
#' @return A list with a plot for each of the transition probability entries
#'   varying with the covariate.
#' @export
#' @import ggplot2
#' @import latex2exp
plot_tpm_entries <- function(num_states, covar_vec, cov_ci,
                             covariate_name = 'Temp - mean(Temp)') {
  plots <- list()
  if (num_states == 3) {
    entries = c(TeX("$\\gamma_{11}$"), TeX("$\\gamma_{21}$"),
                TeX("$\\gamma_{31}$"), TeX("$\\gamma_{12}$"),
                TeX("$\\gamma_{22}$"), TeX("$\\gamma_{32}$"),
                TeX("$\\gamma_{13}$"), TeX("$\\gamma_{23}$"),
                TeX("$\\gamma_{33}$"))
  }
  for (i in 1:(num_states*num_states)) {
    df <- data.frame(Temperature = covar_vec,
                     Estimate = cov_ci$gammas[, i],
                     Lower = cov_ci$lower_gamma[, i],
                     Upper = cov_ci$upper_gamma[, i])
    p <- ggplot(data = df, aes(x = Temperature , y = Estimate)) +
      geom_line() +
      ggtitle(entries[i]) +
      theme(plot.title = ggplot2::element_text(hjust = 0.5),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank())
      #labs(x = covariate_name, y = 'Probability')
    p  <- p + geom_ribbon(data = df, aes(x = Temperature, ymin = Lower,
                                         ymax = Upper), alpha = 0.4)
    plots <- c(plots, list(p))
  }
  plots
}


#' Plot entries of stationary distributions
#'
#' This function plots the stationary distribution entries with their
#' confidence intervals for different covariate values.
#'
#' @param num_states The number of states in the HMM.
#' @param covar_vec The vector of covariate values to plot the entries over.
#' @param entries_list A list.
#' @param covariate_name Label for the x-axis.
#'
#' @return A list with a plot for each of the transition probability entries
#'   varying with the covariate.
#' @export
#' @import ggplot2
#' @import latex2exp
plot_delta_entries <- function(num_states, covar_vec, entries_list,
                               covariate_name = 'Temp - mean(Temp)') {
  plots <- list()
  if (num_states == 3) {
    entries = c(TeX("$\\delta_{1}$"), TeX("$\\delta_{2}$"),
                TeX("$\\delta_{3}$"))
  }
  p <- ggplot(data = df, aes(x = Temperature)) +
    ggtitle(entries[i]) +
    theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
    labs(x = covariate_name, y = 'Probability')
  for (i in 1:(num_states)) {
    df <- data.frame(Temperature = covar_vec,
                     Estimate = entries_list$delta_entries[, i])
      geom_line() +
      ggtitle(entries[i]) +
      theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      labs(x = covariate_name, y = 'Probability')
    #p  <- p + geom_ribbon(data = df, aes(x = Temperature, ymin = Lower,
    #                                     ymax = Upper), alpha = 0.4)
    plots <- c(plots, list(p))
  }
  plots
  #THIS FUNCTION IS WRONG... DO NOT USE
}






