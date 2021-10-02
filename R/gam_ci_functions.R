#' Get indices of parameters in working vector
#'
#' This function finds the indices for each parameter of the gamma HMM within
#' the vector of working parameters outputted by `gam_working_params()`.
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
gam_working_ind <- function(num_states, num_variables, num_subjects,
                             num_covariates, state_dep_dist_pooled = FALSE) {
  alpha_start    <- 1
  alpha_end      <- num_states*num_variables*num_subjects
  theta_end      <- alpha_end + num_states*num_variables*num_subjects
  alpha_len      <- theta_len <- num_subjects*num_states
  if (state_dep_dist_pooled) {
    alpha_end    <- num_states*num_variables
    theta_end    <- alpha_end + num_states*num_variables
    alpha_len    <- theta_len <- num_states
  }
  theta_start <- alpha_end + 1
  beta_start  <- theta_end + 1
  beta_end    <- theta_end + (num_states^2 - num_states)*(num_covariates + 1)
  delta_start <- beta_end + 1
  delta_end   <- delta_start + num_states - 2

  list(alpha_start = alpha_start,
       alpha_end   = alpha_end,
       alpha_len   = alpha_len,
       theta_start = theta_start,
       theta_end   = theta_end,
       theta_len   = theta_len,
       beta_start  = beta_start,
       beta_end    = beta_end,
       delta_start = delta_start,
       delta_end   = delta_end)
}


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
gam_natural_vec <- function(num_states, num_variables, num_subjects,
                            num_covariates, working_params,
                            state_dep_dist_pooled = FALSE) {

  ind <- gam_working_ind(num_states, num_variables, num_subjects,
                         num_covariates, state_dep_dist_pooled = FALSE)

  natural <- numeric()
  natural[ind$alpha_start:ind$alpha_end] <- exp(working_params[ind$alpha_start:
                                                                 ind$alpha_end])
  natural[ind$theta_start:ind$theta_end] <- exp(working_params[ind$theta_start:
                                                                 ind$theta_end])
  natural[ind$beta_start:ind$beta_end]   <- working_params[ind$beta_start:
                                                             ind$beta_end]
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
gam_ci_data <- function(num_states, num_variables, num_subjects,
                        num_covariates, estimate_vec, upper_vec, lower_vec,
                        state_dep_dist_pooled = FALSE) {

  ind <- gam_working_ind(num_states, num_variables, num_subjects,
                         num_covariates, state_dep_dist_pooled = FALSE)

  alpha_estimate <- split_vec(estimate_vec, ind$alpha_start,
                              ind$alpha_end, ind$alpha_len)
  alpha_upper    <- split_vec(upper_vec, ind$alpha_start,
                              ind$alpha_end, ind$alpha_len)
  alpha_lower    <- split_vec(lower_vec, ind$alpha_start,
                              ind$alpha_end, ind$alpha_len)
  theta_estimate <- split_vec(estimate_vec, ind$theta_start,
                              ind$theta_end, ind$theta_len)
  theta_upper    <- split_vec(upper_vec, ind$theta_start,
                              ind$theta_end, ind$theta_len)
  theta_lower    <- split_vec(lower_vec, ind$theta_start,
                              ind$theta_end, ind$theta_len)
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
  alpha          <- list(estimate = alpha_estimate,
                         upper = alpha_upper,
                         lower = alpha_lower)
  theta          <- list(estimate = theta_estimate,
                         upper = theta_upper,
                         lower = theta_lower)

  beta_estimate <- matrix(estimate_vec[ind$beta_start:ind$beta_end],
                          nrow = num_states^2 - num_states)
  beta_upper    <- matrix(upper_vec[ind$beta_start:ind$beta_end],
                          nrow = num_states^2 - num_states)
  beta_lower    <- matrix(lower_vec[ind$beta_start:ind$beta_end],
                          nrow = num_states^2 - num_states)
  beta          <- list(estimate = beta_estimate,
                        upper = beta_upper,
                        lower = beta_lower)

  list(alpha = alpha, theta = theta, beta = beta)
}


#' Compute the confidence intervals
#'
#' This function computes the confidence intervals for the estimated parameters
#' `alpha`, `theta`, and `beta` using the variances outputted from
#' `gam_fit_hmm()` and the Monte Carlo method of estimation.
#'
#' @param hmm A list of parameters that specify the gamma HMM, including
#'   `num_states`, `num_variables`, `num_subjects`,`num_covariates`, `alpha`,
#'   `theta`, `beta`, `delta`.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `alpha` and `theta` should be
#'   treated as equal for all subjects.
#' @param n The number of samples in the Monte Carlo fitting.
#' @param level A number indicating the level of confidence for the desired
#'   interval.
#' @param raw_sample A logical variable indicating wheter to simply output the
#'   `n` sampled `alpha`s and `theta`s from the Monte Carlo estimate.
#'
#' @return Either a list of the sampled `alpha`s and `theta`s or the list of
#'   confidence intervals for each parameter.
#' @export
#' @importFrom stats rgamma quantile
gam_ci <- function(hmm, state_dep_dist_pooled = FALSE, n = 100, level= 0.975,
                    raw_sample = FALSE) {

  num_subjects   <- hmm$num_subjects
  num_states     <- hmm$num_states
  num_variables  <- hmm$num_variables
  num_covariates <- hmm$num_covariates
  working_sd     <- sqrt(diag(hmm$inverse_hessian))
  working_params <- hmm$working_params
  if (num_states != 1) {
    if (state_dep_dist_pooled) {
      d <- num_states - 1
    } else {
      d <- (num_states - 1)*num_subjects
    }
    working_params <- working_params[1:(length(working_params) - d)]
  }
  len_w          <- length(working_sd)
  sample         <- rnorm(len_w, mean = hmm$working_params, sd = working_sd)
  nat            <- gam_natural_vec(num_states, num_variables, num_subjects,
                                    num_covariates, sample,
                                    state_dep_dist_pooled)
  len_n          <- length(nat)
  natural        <- matrix(0, ncol = len_n, nrow = n)
  natural[1, ]   <- nat
  upper          <- numeric()
  lower          <- numeric()

  for (l in 2:n) {
    sample       <- rnorm(len_w, mean = hmm$working_params, sd = working_sd)
    natural[l, ] <- gam_natural_vec(num_states, num_variables, num_subjects,
                                    num_covariates, sample,
                                    state_dep_dist_pooled)
  }
  if (raw_sample) {
    ind <- gam_working_ind(num_states, num_variables, num_subjects,
                           num_covariates, state_dep_dist_pooled = FALSE)
    sample      <- list()
    for (l in 1:n) {
      alpha <- split_vec(natural[l, ], ind$alpha_start,
                         ind$alpha_end, ind$alpha_len)
      theta <- split_vec(natural[l, ], ind$theta_start,
                         ind$theta_end, ind$theta_len)
      for (j in 1:num_variables) {
        alpha[[j]] <- matrix(alpha[[j]], ncol = num_states, byrow = TRUE)
        theta[[j]] <- matrix(theta[[j]], ncol = num_states, byrow = TRUE)
      }
      sample[[l]] <- list(alpha = alpha, theta = theta)
    }
    return(sample)
  }
  for (t in 1:len_n) {
    upper[t] <- quantile(natural[, t], probs = level, na.rm = TRUE)
    lower[t] <- quantile(natural[, t], probs = 1 - level, na.rm = TRUE)
  }
  estimate <- gam_natural_vec(num_states, num_variables, num_subjects,
                              num_covariates, working_params,
                              state_dep_dist_pooled)
  gam_ci_data(num_states, num_variables, num_subjects, num_covariates,
              estimate, upper, lower, state_dep_dist_pooled = FALSE)
}


#' Compute the confidence intervals of fitted gamma distributions
#'
#' This is a helper function for `gam_hist_ci()` which computes the confidence
#' intervals for the fitted gamma state dependent distributions by utilizing
#' `gam_ci()`.
#'
#' @param x The data to be fit with an HMM in the form of a 3D array. The
#'   first index (row) corresponds to time, the second (column) to the
#'   variable number, and the third (matrix number) to the subject number.
#' @param num_states The number of states in the desired HMM.
#' @param num_variables The number of variables in the data.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `alpha` and `theta` should be
#'   treated as equal for all subjects.
#' @param sample A list of `alpha` and `theta` sampled according to `gam_ci()`.
#' @param x_step A value indicating the step length for the range of
#'   observation values.
#' @param n The number of samples in the Monte Carlo fitting.
#' @param level A number indicating the level of confidence for the desired
#'   interval.
#'
#' @return A list containing the upper and lower confidence intervals and the
#'   parameter esimates.
#' @export
#' @importFrom stats dgamma quantile
gam_dist_ci_data <- function(x, num_states, num_variables, num_subjects,
                             sample, state_dep_dist_pooled = FALSE,
                             x_step = 0.2, n = 100, level = 0.975) {
  conf_intervals        <- list()
  for (i in 1:num_subjects) {
    conf_intervals[[i]] <- list()
    s_ind   <- i
    if (state_dep_dist_pooled) {
      s_ind <- 1
    }
    for (j in 1:num_variables) {
      range       <- seq(min(x[, j, i], na.rm = TRUE),
                         max(x[, j, i], na.rm = TRUE), length.out = 100)
      xc          <- length(range)
      density.lst <- list()
      for (k in 1:num_states) {
        densities <- matrix(numeric(xc*n), ncol = xc, nrow = n)
        for (l in 1:n) {
          alpha          <- sample[[l]]$alpha[[j]][s_ind, k]
          theta          <- sample[[l]]$theta[[j]][s_ind, k]
          densities[l, ] <- dgamma(range, shape = alpha, scale = theta)
        }
        density.lst[[k]] <- densities
      }
      upper <- matrix(numeric(xc*num_states), ncol = xc, nrow = num_states)
      lower <- matrix(numeric(xc*num_states), ncol = xc, nrow = num_states)

      for (k in 1:num_states) {
        densities <- density.lst[[k]]
        for (t in 1:xc) {
          upper[k, t] <- quantile(densities[, t], probs = level, na.rm = TRUE)
          lower[k, t] <- quantile(densities[, t], probs = 1 - level,
                                  na.rm = TRUE)
        }
      }
      conf_intervals[[i]][[j]] <- list(range = range,
                                       upper = upper,
                                       lower = lower)
    }
  }
  conf_intervals
}


#' Plot histograms with confidence intervals
#'
#' This function plots the histograms for each subject and variable with the
#' fitted state dependent distributions overlayed and their corresponding
#' confidence intervals.
#'
#' @param x The data to be fit with an HMM in the form of a 3D array. The
#'   first index (row) corresponds to time, the second (column) to the
#'   variable number, and the third (matrix number) to the subject number.
#' @param viterbi A matrix with each column indicating the sequence of states
#'   decoded by `gam_viterbi()` that is supposed to have generated the data of
#'   the subject in the corresponding column.
#' @param num_states The number of states in the desired HMM.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_variables The number of variables in the data.
#' @param hmm A list of parameters that specify the gamma HMM, including
#'   `num_states`, `num_variables`, `num_subjects`, `alpha`, `theta`, `gamma`,
#'   `delta`.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `alpha` and `theta` should be
#'   treated as equal for all subjects.
#' @param width The width of the histogram bins.
#' @param n The number of samples in the Monte Carlo fitting.
#' @param level A number indicating the level of confidence for the desired
#'   interval.
#' @param x_step A value indicating the step length for the range of
#'   observation values.
#'
#' @return Histograms of the data with overlayed distributions and confidence
#'   intervals.
#' @export
#' @import RColorBrewer
#' @importFrom stats dgamma
#' @importFrom ggplot2 ggplot geom_histogram aes theme_bw geom_ribbon geom_line
#'   ggtitle theme labs
gam_hist_ci <- function(x, viterbi, num_states, num_subjects, num_variables,
                        hmm, state_dep_dist_pooled = FALSE,
                        width = 1, n = 100, level = 0.975, x_step = 0.2) {

  sample         <- gam_ci(hmm, state_dep_dist_pooled, n,
                           level, raw_sample = TRUE)
  conf_intervals <- gam_dist_ci_data(x, num_states, num_variables,
                                     num_subjects, sample,
                                     state_dep_dist_pooled, x_step, n, level)
  n     <- nrow(x)
  Var   <- c("Variable 1", "Variable 2", "Variable 3", "Variable 4")
  Sub   <- c("Subject 1", "Subject 2", "Subject 3", "Subject 4")
  plots <- list()
  for (i in 1:num_subjects) {
    subvar_data <- data.frame('State' = as.factor(viterbi[, i]))
    s_ind   <- i
    if (state_dep_dist_pooled) {
      s_ind <- 1
    }
    for (j in 1:num_variables) {
      subvar_data$Observation <- x[, j, i]
      h <- ggplot2::ggplot() +
        ggplot2::geom_histogram(data = subvar_data,
                                aes(x = Observation),
                                binwidth = width,
                                colour = "grey",
                                fill = "white") +
        ggplot2::scale_color_brewer(palette = "Set1") +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(Sub[i]) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(x = Var[j], y = '')

      xfit <- seq(min(subvar_data$Observation, na.rm = TRUE),
                  max(subvar_data$Observation, na.rm = TRUE),
                  by = x_step)
      marginal <- numeric(length(xfit))
      for (k in 1:num_states) {
        yfit     <- dgamma(xfit, shape = hmm$alpha[[j]][s_ind, k],
                           scale = hmm$theta[[j]][s_ind, k])
        yfit     <- yfit * sum(subvar_data$State == k) * width
        df       <- data.frame('xfit' = xfit, 'yfit' = yfit,
                               col = as.factor(rep(k, length(xfit))))
        h        <- h + ggplot2::geom_line(data = df,
                                           aes(xfit, yfit, colour = col),
                                           lwd = 0.7)
        marginal <- marginal + yfit
      }
      h  <- h + labs(color = "State")
      df <- data.frame('xfit' = xfit, 'yfit' = marginal)
      h  <- h + geom_line(data = df, aes(xfit, yfit), col = "black", lwd = 0.7)

      for (k in 1:num_states){
        upper <- conf_intervals[[i]][[j]]$upper[k, ]*
          sum(subvar_data$State == k)*width
        lower <- conf_intervals[[i]][[j]]$lower[k, ]*
          sum(subvar_data$State == k)*width
        df <- data.frame('x' = conf_intervals[[i]][[j]]$range,
                            'upper' = upper, 'lower' = lower)
        h <- h + ggplot2::geom_ribbon(data = df,
                                      aes(x = x, ymin = lower, ymax = upper),
                                      alpha = 0.4,
                                      fill = c(RColorBrewer::brewer.pal(n = 8,
                                                          name = "Set1"))[k])
      }
      plots <- c(plots, list(h))
    }
  }
  plots
  # ggarrange(plotlist = plots, common.legend = TRUE, legend = "bottom",
  #           labels = c("Subject 1",
  #                      "Subject 1",
  #                      "Subject 2",
  #                      "Subject 2"))
}
