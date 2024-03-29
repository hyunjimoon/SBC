% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\name{generate_datasets}
\alias{generate_datasets}
\title{Generate datasets.}
\usage{
generate_datasets(generator, n_sims, n_datasets = NULL)
}
\arguments{
\item{generator}{a generator object - build e.g. via \code{SBC_generator_function} or
\code{SBC_generator_brms}.}

\item{n_sims}{the number of simulated datasets to use}

\item{n_datasets}{DEPRECATED, use \code{n_sims} instead.}
}
\value{
object of class \code{SBC_datasets}
TODO: seed
}
\description{
Generate datasets.
}
