#' Compute transition probability matrices for plotting
#'
#' This function computes the transition probability matrix entries with their
#' confidence intervals for different covariate values in order to plot how
#' they change.
#'
#' @param num_states The number of states in the desired HMM.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_covariates The number of covariates in the data that the
#'   transition probability matrix depends on.
#' @param num_time A value indicating the length of the data.
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
tpm_entries <- function(num_states, num_subjects, num_covariates, num_time,
                        design, conf_intervals) {
  estimate <- fit_tpm(num_states, num_subjects, num_covariates, num_time,
                      conf_intervals$beta$estimate, design)
  estimate_entries <- matrix(unlist(estimate), ncol = num_states*num_states,
                             byrow = TRUE)
  #lower <- fit_tpm(num_states, num_subjects, num_covariates, num_time,
  #                 conf_intervals$beta$lower, design)
  #lower_entries <- matrix(unlist(lower), ncol = num_states*num_states,
  #                        byrow = TRUE)
  #upper <- fit_tpm(num_states, num_subjects, num_covariates, num_time,
  #                 conf_intervals$beta$upper, design)
  #upper_entries <- matrix(unlist(upper), ncol = num_states*num_states,
  #                        byrow = TRUE)
  delta_list <- list()
  for (i in num_subjects) {
    delta_list[[i]] <- matrix(0, ncol = num_states, nrow = num_time)
    for (t in 1:num_time){
      eigs <- eigen(t(estimate[[i]][, , t]))
      values <- numeric()
      for (value in eigs$value) {
        if (zapsmall(Im(value)) == 0) {
          values <- c(values, zapsmall(Re(value)))
        } else {
          values <- c(values, zapsmall(value))
        }
      }
      ind <- which(values == 1)
      d <- eigs$vectors[, ind]/sum(eigs$vectors[, ind])
      delta_list[[i]][t, ] <- eigs$vectors[, ind]/sum(eigs$vectors[, ind])
      #return(list(d, d %*% estimate[[i]][, , t]))
    }
  }
  delta_entries = matrix(unlist(delta_list), ncol = num_states, byrow = FALSE)
  list(estimate_entries = estimate_entries,
       delta_entries = delta_entries)
       #lower_entries = lower_entries,
       #upper_entries = upper_entries)
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
plot_tpm_entries <- function(num_states, covar_vec, tpm_entry_list,
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
                     Estimate = tpm_entry_list$estimate_entries[, i])
                     #Lower = tpm_entry_list$lower_entries[, i],
                     #Upper = tpm_entry_list$upper_entries[, i])
    p <- ggplot(data = df, aes(x = Temperature , y = Estimate)) +
      geom_line() +
      ggtitle(entries[i]) +
      theme(plot.title = ggplot2::element_text(hjust = 0.5),
            axis.title.x = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank())
      #labs(x = covariate_name, y = 'Probability')
    #p  <- p + geom_ribbon(data = df, aes(x = Temperature, ymin = Lower,
    #                                     ymax = Upper), alpha = 0.4)
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






