% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deconvolution.R
\name{.estimate_contingency_table}
\alias{.estimate_contingency_table}
\title{Fast deconvolution of read count contingency table. Used as basis of
the exported fast_estimate_* functions}
\usage{
.estimate_contingency_table(
  total_readdepth,
  alt_readdepth,
  logr,
  host_total_readdepth,
  host_alt_readdepth,
  purity,
  ploidy,
  host_ploidy = 2
)
}
\arguments{
\item{total_readdepth}{Tumour sample total read depth}

\item{alt_readdepth}{Tumour sample alt read depth}

\item{logr}{Tumour sample logR}

\item{host_total_readdepth}{Host sample total read depth}

\item{host_alt_readdepth}{Host sample alt read depth}

\item{purity}{Estimated purity of tumour sample
       |  Ref  |  Alt  |
-------|-------|-------|-------
 Host  |  (a)  |  (b)  |  (c)
Tumour |  (d)  |  (e)  |  (f)
-------|-------|-------|-------
       |   R   |   A   |   T
   |  Ref  |  Alt  |
---|-------|-------|-------
 H |  K-L  |   L   |   K
 T |T-K-A+L|  A-L  |  T-K
---|-------|-------|-------
   |   R   |   A   |   T}
}
\description{
Fast deconvolution of read count contingency table. Used as basis of
the exported fast_estimate_* functions
}
