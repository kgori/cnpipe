% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ascat_functions.R
\name{winsorize}
\alias{winsorize}
\title{Winsorizes an input vector. Winsorizing means to clip outlier values to a level closer to the central tendency of the distribution.}
\usage{
winsorize(x, tau = 2.5, k = 25)
}
\arguments{
\item{tau}{Multiplier of the local standard deviation: smaller values result in a stronger clipping effect}

\item{k}{The local window size}
}
\description{
Winsorizes an input vector. Winsorizing means to clip outlier values to a level closer to the central tendency of the distribution.
}
