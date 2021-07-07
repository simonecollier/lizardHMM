#' Plot timeseries of observations
#'
#' @inheritParams norm_fit_hmm
#' @param states A vector containing the sequence of states that generated the
#'   data `x`.
#' @param variable_names A vector containing the names of the variables in the
#'   data `x`.
#' @param subject_names A vector containing the names of the subjects generating
#'   the data `x`.
#' @param length A number indicating the desired length of the timeseries plot.
#'
#' @return A grid of plots of the time series for each subject and variable.
#' @export
#' @importFrom ggplot2 ggplot aes ggtitle theme labs geom_point geom_line

timeseries_plot <- function(x, states, num_subjects, num_variables,
                            variable_names = c("Var 1", "Var 2", "Var 3"),
                            subject_names = c("Subject 1", "Subject 2",
                                              "Subject 3", "Subject 4"),
                            length = 300) {
  n      <- nrow(x)
  plots  <- list()
  for (i in 1:num_subjects) {
    data       <- data.frame('State' = as.factor(states[1:length, i]))
    data$Time  <- 1:length
    for (j in 1:num_variables) {
      data$Observation <- x[1:length, j, i]
      p <- ggplot(data, ggplot2::aes(x = Time, y = Observation)) +
        ggplot2::theme_light() +
        ggtitle(subject_names[i]) +
        theme(axis.title.x = ggplot2::element_blank(),
              plot.title = ggplot2::element_text(hjust = 0.5)) +
        labs(x = '', y = variable_names[j]) +
        geom_point(ggplot2::aes(color = State)) +
        geom_line(colour = 'grey', alpha = 0.8, lwd = 0.4)
      plots <- c(plots, list(p))
    }
  }
  plots
  #egg::ggarrange(plotlist = plots, common.legend = TRUE, legend = "bottom")
}
