#' Plot normal pseudo-residuals
#'
#' This function produces a scatter plot of the normal pseudo-residuals in time
#' for each subject.
#'
#' @param forecast_psr A list containing a vector of pseudo-residuals for each
#'   subject.
#' @param num_subjects The number of subjects/trials that generated the data.
#'
#' @return A plot of the normal-pseudo residuals for each subject.
#' @export
#' @importFrom ggplot2 ggplot aes theme theme_bw ylab geom_point geom_hline
psr_plot <- function(forecast_psr, num_subjects) {
  plots <- list()
  for (i in num_subjects) {
    data <- data.frame('Time' = 1:length(forecast_psr[[i]]),
                       'pseudo_residual' = forecast_psr[[i]])
    p <- ggplot(data) +
      geom_point(aes(x = Time, y = pseudo_residual), size = 0.04, alpha = 0.6) +
      theme_bw() +
      theme(panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()) +
      geom_hline(yintercept = 0, linetype = 2, alpha = 0.9) +
      ylab("Normal Pseudo Residuals")
    plots <- c(plots, list(p))
  }
  # if (labs==FALSE){
  #   p <- p + theme(panel.grid.major = element_blank(),
  #                  panel.grid.minor = element_blank(),
  #                  axis.title.x = element_blank(),
  #                  axis.title.y = element_blank())
  # }
  plots
}


#' Plot normal pseudo-residuals
#'
#' This function produces a histogram of the normal pseudo-residuals for each
#' subject.
#'
#' @param forecast_psr A list containing a vector of pseudo-residuals for each
#'   subject.
#' @param num_subjects The number of subjects/trials that generated the data.
#'
#' @return A histogram of the normal-pseudo residuals for each subject.
#' @export
#' @importFrom ggplot2 ggplot aes theme theme_bw xlab xlim stat_function
#'   geom_histogram
psr_hist <- function(forecast_psr, num_subjects) {
  plots <- list()
  xaxis <- seq(-4, 4, 0.01)
  # norm_data <- pnorm(seq(-4, 4, 0.01))^4
  # norm_data <- qnorm(norm_data)
  for (i in num_subjects) {
    data <- data.frame('Time' = 1:length(forecast_psr[[i]]),
                       'pseudo_residual' = forecast_psr[[i]])
    p <- ggplot(data, aes(pseudo_residual)) +
      geom_histogram(aes(y = ..density..), binwidth = 0.5,
                     colour="white", fill="grey") +
      stat_function(fun = dnorm, colour = "black") +
      #geom_line(x = xaxis, y = norm_data, colour = "black") +
      theme_bw() +
      theme(panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank()) +
      xlab('normal pseudo-residuals') +
      xlim(c(-4, 4))
    plots <- c(plots, list(p))
  }
  # if (labs==FALSE){
  #   p <- p + theme(panel.grid.major = element_blank(),
  #                  panel.grid.minor = element_blank(),
  #                  axis.title.x = element_blank(),
  #                  axis.title.y = element_blank())
  # }
  plots
}


#' Quantile-quantile plot of normal pseudo-residuals
#'
#' This function produces a quantile-quantile plot of the normal
#' pseudo-residuals for each subject.
#'
#' @param forecast_psr A list containing a vector of pseudo-residuals for each
#'   subject.
#' @param num_subjects The number of subjects/trials that generated the data.
#'
#' @return A qq plot of the normal-pseudo residuals for each subject.
#' @export
#' @importFrom ggplot2 ggplot aes theme theme_bw stat_qq stat_qq_line
psr_qq <- function(forecast_psr, num_subjects) {
  plots <- list()
  for (i in num_subjects) {
    data <- data.frame('Time' = 1:length(forecast_psr[[i]]),
                       'pseudo_residual' = forecast_psr[[i]])
    p <- ggplot(data, aes(sample = pseudo_residual, alpha = 0.5)) +
      stat_qq(size = 0.1) +
      stat_qq_line() +
      theme_bw() +
      theme(panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            legend.position = "none")
    plots <- c(plots, list(p))
  }
  plots
}


#' ACF plot of normal pseudo-residuals
#'
#' This function produces a ACF plot of the normal pseudo-residuals for each
#' subject.
#'
#' @param forecast_psr A list containing a vector of pseudo-residuals for each
#'   subject.
#' @param num_subjects The number of subjects/trials that generated the data.
#'
#' @return An acf plot of the normal-pseudo residuals for each subject.
#' @export
#' @importFrom ggplot2 ggplot aes theme theme_bw
psr_acf <- function(forecast_psr, num_subjects) {
  plots <- list()
  for (i in num_subjects) {
    p <- forecast::ggAcf(x = forecast_psr[[i]]) +
      theme_bw() +
      theme(panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            plot.title = ggplot2::element_blank())
    plots <- c(plots, list(p))
  }
  plots
}
