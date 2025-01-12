% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binomial_mixture_model.R
\name{binomial_MM}
\alias{binomial_MM}
\title{Fits a binomial mixture model to data in y. The number of components
is inferred from the length of the theta starting parameters vector.}
\usage{
binomial_MM(y, theta, lambda, tol = 0.001, max_iter = 1000)
}
\arguments{
\item{y}{Data to fit, as matrix of successes and failures in columns
(in DNA variant context, a 'success' is an ALT read, so they go in
column 1)}

\item{theta}{Vector of starting binomial probability parameters,
one for each mixture component}

\item{lambda}{Vector of mixture weights, one for each component}

\item{tol}{Optimisation will stop when the bmm_likelihood improves by
less than this tolerance}

\item{max_iter}{Optimisation stops when this number of iterations is
reached}
}
\description{
Fits a binomial mixture model to data in y. The number of components
is inferred from the length of the theta starting parameters vector.
}
