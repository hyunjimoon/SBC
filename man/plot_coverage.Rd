% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_coverage}
\alias{plot_coverage}
\alias{plot_coverage_diff}
\title{Plot the observed coverage and its uncertainty.}
\usage{
plot_coverage(
  x,
  variables = NULL,
  prob = 0.95,
  interval_type = "central",
  parameters = NULL,
  max_points = NULL
)

plot_coverage_diff(
  x,
  variables = NULL,
  prob = 0.95,
  interval_type = "central",
  parameters = NULL,
  max_points = NULL
)
}
\arguments{
\item{x}{object containing results (a data.frame or \code{\link[=SBC_results]{SBC_results()}} object).}

\item{variables}{variables to show in the plot or \code{NULL} to show all}

\item{prob}{the with of the uncertainty interval to be shown}

\item{parameters}{DEPRECATED. Use \code{variables} instead.}

\item{max_points}{maximum number of points where to evaluate the coverage.
If set to \code{NULL}, coverage is evaluated across the whole range of ranks.
Setting to some smaller number may reduce memory footprint and increase speed.}
}
\value{
a ggplot2 plot object
}
\description{
\code{plot_coverage} will plot the observed coverage,
while \code{plot_coverage_diff} will show the difference between observed
and expected coverage.
Please refer to \code{\link[=empirical_coverage]{empirical_coverage()}} for details on computation
and limitations of this plot as well as details on the arguments.
See \code{vignette("rank_visualizations")} for
more details.
}
\seealso{
empirical_coverage
}
