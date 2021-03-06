% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covariate_analysis.R
\name{plot_delta_entries}
\alias{plot_delta_entries}
\title{Plot entries of stationary distributions}
\usage{
plot_delta_entries(
  num_states,
  covar_vec,
  entries_list,
  covariate_name = "Temp - mean(Temp)"
)
}
\arguments{
\item{num_states}{The number of states in the HMM.}

\item{covar_vec}{The vector of covariate values to plot the entries over.}

\item{entries_list}{A list.}

\item{covariate_name}{Label for the x-axis.}
}
\value{
A list with a plot for each of the transition probability entries
varying with the covariate.
}
\description{
This function plots the stationary distribution entries with their
confidence intervals for different covariate values.
}
