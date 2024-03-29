% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metric.R
\name{cjs_dist}
\alias{cjs_dist}
\title{Cumulative Jensen-Shannon divergence}
\usage{
cjs_dist(
  x,
  y,
  x_weights = rep(1/length(x), length(x)),
  y_weights = rep(1/length(y), length(y)),
  ...
)
}
\arguments{
\item{x}{numeric vector of draws from first distribution}

\item{y}{numeric vector of draws from second distribution}

\item{x_weights}{numeric vector of weights of first distribution}

\item{y_weights}{numeric vector of weights of second distribution}

\item{...}{unused}
}
\value{
distance value based on CJS computation.
}
\description{
Computes the cumulative Jensen-Shannon distance between two
samples.
}
\details{
The Cumulative Jensen-Shannon distance is a symmetric metric based
on the cumulative Jensen-Shannon divergence. The divergence CJS(P || Q) between
two cumulative distribution functions P and Q is defined as:

\deqn{CJS(P || Q) = \sum P(x) \log \frac{P(x)}{0.5 (P(x) + Q(x))} + \frac{1}{2 \ln 2} \sum (Q(x) - P(x))}

The symmetric metric is defined as:

\deqn{CJS_{dist}(P || Q) = \sqrt{CJS(P || Q) + CJS(Q || P)}}

This has an upper bound of \eqn{\sqrt \sum (P(x) + Q(x))}
}
\references{
Nguyen H-V., Vreeken J. (2015).  Non-parametric
Jensen-Shannon Divergence.  In: Appice A., Rodrigues P., Santos
Costa V., Gama J., Jorge A., Soares C. (eds) Machine Learning
and Knowledge Discovery in Databases.  ECML PKDD 2015. Lecture
Notes in Computer Science, vol 9285.  Springer, Cham.
\code{doi:10.1007/978-3-319-23525-7_11}
}
