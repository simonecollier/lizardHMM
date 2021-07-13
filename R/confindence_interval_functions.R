#' Convert working parameters to a vector of natural parameters
#'
#' @inheritParams norm_natural_params
#'
#' @return A vector of the natural parameters.
#' @export
#'
#' @examples
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
#' @inheritParams natural_params
#' @param estimate_vec A vector containing the estimated natural parameters.
#' @param upper_vec A vector containing the upper confidence interval of the
#'   estimated natural parameters.
#' @param lower_vec A vector containing the lower confidence interval of the
#'   estimated natural parameters.
#'
#' @return A list containing the upper and lower confidence interval for each
#'   estimated parameter.
#' @export
#'
#' @examples
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
#' This function computes the standard errors of the fitted HMM parameters
#' using the inverse hessian computed from `norm_fit_hmm()` and the Monte Carlo
#' method.
#'
#' @inheritParams norm_generate_sample
#' @param n The number of samples in the Monte Carlo fitting.
#' @param level A number indicating the level of confidence for the desired
#'   interval.
#' @param raw_sample A logical vaiable indicating whether to output just
#'   the raw samples without the computed confidence intervals.
#'
#' @return Either a matrix of the samples natural parameters or the list of
#'   confidence interals for each parameter.
#' @export
#' @importFrom stats rnorm quantile
#'
#' @examples
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
    return(natural)
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


#' #' Compute the confidence intervals of fitted normal distributions
#'
#' This function computes the confidence intervals for the fitted normal state
#' dependent distributions using the Monte Carlo approach.
#'
#' @inheritParams norm_loglikelihood
#' @inheritParams norm_ci
#' @param sample_params A matrix of the natural parameters.
#' @param x_step A value indicating the step length for the range of
#'   observation values.
#'
#' @return
#' @export
#' @importFrom stats dnorm quantile
norm_dist_ci_data <- function(x, num_states, num_variables, num_subjects,
                              sample_params, state_dep_dist_pooled = FALSE,
                              x_step = 0.2, n = 100, level = 0.975) {
  conf_intervals <- list()
  ns             <- num_subjects
  if (state_dep_dist_pooled) {
    ns <- 1
  }
  mu_ind         <- 1
  sigma_ind      <- ns*num_variables*num_states
  for (j in 1:num_variables) {
    for (i in 1:ns) {
      conf_intervals[[i]] <- list()
      range       <- seq(min(x[, j, i]), max(x[, j, i]), length.out = 100)
      xc          <- length(range)
      density.lst <- list()
      for (k in 1:num_states) {
        densities <- matrix(numeric(xc*n), ncol = xc, nrow = n)
        mu_ind    <- mu_ind + 1
        for (l in 1:n) {
          mu             <- sample_params[l, mu_ind]
          sigma          <- sample_params[l, (sigma_ind + mu_ind)]
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


#' Title
#'
#' @param x
#' @param viterbi
#' @param conf_intervals
#' @param num_states
#' @param num_subjects
#' @param num_variables
#' @param hmm
#' @param state_dep_dist_pooled
#' @param width
#' @param n
#' @param level
#' @param x_step
#'
#' @return
#' @export
#' @importFrom stats dnorm
#' @importFrom ggplot2 ggplot geom_histogram aes theme_bw geom_ribbon geom_line
#'   ggtitle theme labs
#'
#' @examples
norm_hist_ci <- function(x, viterbi, num_states, num_subjects, num_variables,
                         hmm, state_dep_dist_pooled = FALSE,
                         width = 1, n = 100, level = 0.975, x_step = 0.2) {

  natural        <- norm_ci(hmm, state_dep_dist_pooled, n,
                            level, raw_sample = TRUE)
  conf_intervals <- norm_dist_ci_data(x, num_states, num_variables,
                                      num_subjects, natural,
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
      h <- ggplot() +
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
