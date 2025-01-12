% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvolution.R
\name{prob_read_came_from_host}
\alias{prob_read_came_from_host}
\title{Probability that a read came from the host. See `prob_read_came_from_tumour`}
\usage{
prob_read_came_from_host(logr, purity, ploidy, host_ploidy)
}
\arguments{
\item{logr}{Tumour sample logR}

\item{purity}{Tumour sample purity}

\item{ploidy}{Tumour ploidy (ploidy of pure tumour)}

\item{host_ploidy}{Host ploidy estimate (almost always will be 2)}

\item{host_copynumber}{Host copy number estimate (usually 2)}
}
\description{
Probability that a read came from the host. See `prob_read_came_from_tumour`
}
