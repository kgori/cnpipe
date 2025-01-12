% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ascat_functions.R
\name{calc_total_copynumber}
\alias{calc_total_copynumber}
\title{Estimates the tumour copy number from the logR value.}
\usage{
calc_total_copynumber(
  logr,
  purity,
  ploidy,
  host_copy_number = 2,
  host_ploidy = 2
)
}
\arguments{
\item{logr}{- log2 of the ratio between tumour and host read depths. The read depths should be normalised prior
to computing the logr, so that any difference in sequencing depth between the host and tumour samples doesn't
influence the estimate of copy number.}

\item{purity}{- Fraction of the tumour sample that is derived from tumour cells.}

\item{ploidy}{- Overall tumour ploidy. The average read depth over the genome corresponds to this copy number.}

\item{host_copy_number}{- The estimated copy number in the host for the genomic region being considered.}

\item{host_ploidy}{- Overall host ploidy. (default = 2)}
}
\description{
Estimates the tumour copy number from the logR value.
}
