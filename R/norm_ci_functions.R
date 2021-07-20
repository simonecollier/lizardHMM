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
#' @param working_params A vector of the working normal parameters for the
#'   HMM.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
#'   treated as equal for all subjects.
#'
#' @return A vector of the natural parameters.
#' @export
norm_natural_vec <- function(num_states, num_variables, num_subjects,
                             num_covariates, working_params,
                             state_dep_dist_pooled = FALSE) {
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

  natural                        <- numeric()
  natural[mu_start:mu_end]       <- working_params[mu_start:mu_end]
  natural[sigma_start:sigma_end] <- exp(working_params[sigma_start:sigma_end])
  natural[beta_start:beta_end]   <- working_params[beta_start:beta_end]

  natural
}


#' Reformat confidence interval data
#'
#' This is a helper function for `norm_ci()` which reformats the output so that
#' it is more easily interpreted.
#'
#' @param num_states The number of states in the desired HMM.
#' @param num_variables The number of variables in the data.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_covariates The number of covariates in the data that the
#'   transition probability depends on.
#' @param estimate_vec A vector containing the estimated natural parameters
#'   of the normal HMM in the format outputte by `norm_natural_vec()`.
#' @param upper_vec A vector containing the upper confidence interval of the
#'   estimated natural parameters of the normal HMM in the format outputte by
#'   `norm_natural_vec()`.
#' @param lower_vec A vector containing the lower confidence interval of the
#'   estimated natural parameters of the normal HMM in the format outputte by
#'   `norm_natural_vec()`.
#'@param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
#'   treated as equal for all subjects.
#'
#' @return A list containing the upper and lower confidence interval for each
#'   estimated parameter.
#' @export
norm_ci_data <- function(num_states, num_variables, num_subjects,
                         num_covariates, estimate_vec, upper_vec, lower_vec,
                         state_dep_dist_pooled = FALSE) {
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

  mu_estimate    <- split_vec(estimate_vec, mu_start, mu_end, mu_len)
  mu_upper       <- split_vec(upper_vec, mu_start, mu_end, mu_len)
  mu_lower       <- split_vec(lower_vec, mu_start, mu_end, mu_len)
  sigma_estimate <- split_vec(estimate_vec, sigma_start, sigma_end, sigma_len)
  sigma_upper    <- split_vec(upper_vec, sigma_start, sigma_end, sigma_len)
  sigma_lower    <- split_vec(lower_vec, sigma_start, sigma_end, sigma_len)
  for (j in 1:num_variables) {
    mu_estimate[[j]]    <- matrix(mu_estimate[[j]], ncol = num_states,
                                  byrow = TRUE)
    mu_upper[[j]]       <- matrix(mu_upper[[j]], ncol = num_states,
                                  byrow = TRUE)
    mu_lower[[j]]       <- matrix(mu_lower[[j]], ncol = num_states,
                                  byrow = TRUE)
    sigma_estimate[[j]] <- matrix(sigma_estimate[[j]], ncol = num_states,
                                  byrow = TRUE)
    sigma_upper[[j]]    <- matrix(sigma_upper[[j]], ncol = num_states,
                                  byrow = TRUE)
    sigma_lower[[j]]    <- matrix(sigma_lower[[j]], ncol = num_states,
                                  byrow = TRUE)
  }
  mu             <- list(estimate = mu_estimate,
                         upper = mu_upper,
                         lower = mu_lower)
  sigma          <- list(estimate = sigma_estimate,
                         upper = sigma_upper,
                         lower = sigma_lower)

  beta_estimate <- matrix(estimate_vec[beta_start:beta_end],
                          nrow = num_states^2 - num_states)
  beta_upper    <- matrix(upper_vec[beta_start:beta_end],
                          nrow = num_states^2 - num_states)
  beta_lower    <- matrix(lower_vec[beta_start:beta_end],
                          nrow = num_states^2 - num_states)
  beta          <- list(estimate = beta_estimate,
                        upper = beta_upper,
                        lower = beta_lower)

  list(mu = mu, sigma = sigma, beta = beta)
}


#' Compute the confidence intervals
#'
#' This function computes the confidence intervals for the estimated parameters
#' `mu`, `sigma`, and `beta` using the variances outputted from `norm_fit_hmm()`
#' and the Monte Carlo method of estimation.
#'
#' @param hmm A list of parameters that specify the normal HMM, including
#'   `num_states`, `num_variables`, `num_subjects`,`num_covariates`, `mu`,
#'   `sigma`, `beta`, `delta`.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
#'   treated as equal for all subjects.
#' @param n The number of samples in the Monte Carlo fitting.
#' @param level A number indicating the level of confidence for the desired
#'   interval.
#' @param raw_sample A logical variable indicating wheter to simply output the
#'   `n` sampled `mu`s and `sigma`s from the Monte Carlo estimate.
#'
#' @return Either a list of the sampled `mu`s and `sigma`s or the list of
#'   confidence intervals for each parameter.
#' @export
#' @importFrom stats rnorm quantile
norm_ci <- function(hmm, state_dep_dist_pooled = FALSE, n = 100, level= 0.975,
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
  nat            <- norm_natural_vec(num_states, num_variables, num_subjects,
                                     num_covariates, sample,
                                     state_dep_dist_pooled)
  len_n          <- length(nat)
  natural        <- matrix(0, ncol = len_n, nrow = n)
  natural[1, ]   <- nat
  upper          <- numeric()
  lower          <- numeric()

  for (l in 2:n) {
    sample       <- rnorm(len_w, mean = hmm$working_params, sd = working_sd)
    natural[l, ] <- norm_natural_vec(num_states, num_variables, num_subjects,
                                     num_covariates, sample,
                                     state_dep_dist_pooled)
  }
  if (raw_sample) {
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
    sample      <- list()
    for (l in 1:n) {
      mu    <- split_vec(natural[l, ], mu_start, mu_end, mu_len)
      sigma <- split_vec(natural[l, ], sigma_start, sigma_end, sigma_len)
      for (j in 1:num_variables) {
        mu[[j]]    <- matrix(mu[[j]], ncol = num_states, byrow = TRUE)
        sigma[[j]] <- matrix(sigma[[j]], ncol = num_states, byrow = TRUE)
      }
      sample[[l]] <- list(mu = mu, sigma = sigma)
    }
    return(sample)
  }
  for (t in 1:len_n) {
    upper[t] <- quantile(natural[, t], probs = level, na.rm = TRUE)
    lower[t] <- quantile(natural[, t], probs = 1 - level, na.rm = TRUE)
  }
  estimate <- norm_natural_vec(num_states, num_variables, num_subjects,
                               num_covariates, working_params,
                               state_dep_dist_pooled)
  norm_ci_data(num_states, num_variables, num_subjects, num_covariates,
               estimate, upper, lower, state_dep_dist_pooled = FALSE)
}


#' Compute the confidence intervals of fitted normal distributions
#'
#' This is a helper function for `norm_hist_ci()` which computes the confidence
#' intervals for the fitted normal state dependent distributions by utilizing
#' `norm_ci()`.
#'
#' @param x The data to be fit with an HMM in the form of a 3D array. The
#'   first index (row) corresponds to time, the second (column) to the
#'   variable number, and the third (matrix number) to the subject number.
#' @param num_states The number of states in the desired HMM.
#' @param num_variables The number of variables in the data.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
#'   treated as equal for all subjects.
#' @param sample A list of `mu` and `sigma` sampled according to `norm_ci()`.
#' @param x_step A value indicating the step length for the range of
#'   observation values.
#' @param n The number of samples in the Monte Carlo fitting.
#' @param level A number indicating the level of confidence for the desired
#'   interval.
#'
#' @return A list containing the upper and lower confidence intervals and the
#'   parameter esimates.
#' @export
#' @importFrom stats dnorm quantile
norm_dist_ci_data <- function(x, num_states, num_variables, num_subjects,
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
      range       <- seq(min(x[, j, i]), max(x[, j, i]), length.out = 100)
      xc          <- length(range)
      density.lst <- list()
      for (k in 1:num_states) {
        densities <- matrix(numeric(xc*n), ncol = xc, nrow = n)
        for (l in 1:n) {
          mu             <- sample[[l]]$mu[[j]][s_ind, k]
          sigma          <- sample[[l]]$sigma[[j]][s_ind, k]
          densities[l, ] <- dnorm(range, mu, sigma)
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
#'   decoded by `norm_viterbi()` that is supposed to have generated the data of
#'   the subject in the corresponding column.
#' @param num_states The number of states in the desired HMM.
#' @param num_subjects The number of subjects/trials that generated the data.
#' @param num_variables The number of variables in the data.
#' @param hmm A list of parameters that specify the normal HMM, including
#'   `num_states`, `num_variables`, `num_subjects`, `mu`, `sigma`, `gamma`,
#'   `delta`.
#' @param state_dep_dist_pooled A logical variable indiacting whether the
#'   state dependent distribution parameters `mu` and `sigma` should be
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
#' @importFrom stats dnorm
#' @importFrom ggplot2 ggplot geom_histogram aes theme_bw geom_ribbon geom_line
#'   ggtitle theme labs
norm_hist_ci <- function(x, viterbi, num_states, num_subjects, num_variables,
                         hmm, state_dep_dist_pooled = FALSE,
                         width = 1, n = 100, level = 0.975, x_step = 0.2) {

  sample        <- norm_ci(hmm, state_dep_dist_pooled, n,
                            level, raw_sample = TRUE)
  conf_intervals <- norm_dist_ci_data(x, num_states, num_variables,
                                      num_subjects, sample,
                                      state_dep_dist_pooled, x_step, n, level)
  n      <- nrow(x)
  Var    <- c("Variable 1", "Variable 2", "Variable 3", "Variable 4")
  Sub    <- c("Subject 1", "Subject 2", "Subject 3", "Subject 4")
  plots  <- list()
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
                                colour = "cornsilk4",
                                fill = "white") +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(Sub[i]) +
        ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
                       panel.grid.minor = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(hjust = 0.5)) +
        ggplot2::labs(x = Var[j], y = '')

      xfit <- seq(min(subvar_data$Observation), max(subvar_data$Observation),
                  by = x_step)
      marginal <- numeric(length(xfit))
      for (k in 1:num_states) {
        yfit     <- dnorm(xfit, hmm$mu[[j]][s_ind, k], hmm$sigma[[j]][s_ind, k])
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
      h  <- h + geom_line(data = df, aes(xfit, yfit), col="black", lwd=0.7)

      for (k in 1:num_states){
        upper <- conf_intervals[[i]][[j]]$upper[k, ]*
          sum(subvar_data$State == k)*width
        lower <- conf_intervals[[i]][[j]]$lower[k, ]*
          sum(subvar_data$State == k)*width
        df <- data.frame('x' = conf_intervals[[i]][[j]]$range,
                         'upper' = upper, 'lower' = lower)
        h <- h + ggplot2::geom_ribbon(data = df,
                                      aes(x = x, ymin = lower, ymax = upper),
                                      fill = (k + 1), alpha = 0.4)
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
