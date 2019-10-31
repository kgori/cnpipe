#' Calculate the log odds ratio (logOR) of tumour and reference allele counts,
#' as an alternative to VAF. Formulation is from Facets package.
#' @param tumour_alt - count of tumour alt reads at a site
#' @param tumour_ref - count of tumour ref reads at a site
#' @param host_alt - count of host alt reads at a site
#' @param host_ref - count of host ref reads at a site
log_odds_ratio <- function(tumour_alt, tumour_ref, host_alt, host_ref) {
    log2(tumour_alt + 0.5) - log2(tumour_ref + 0.5) - log2(host_alt + 0.5) + log2(host_ref + 0.5)
}

#' Calculate the variance of the log odds ratio (logOR) of tumour and reference allele counts.
#' Formulation is from Facets package.
#' @param tumour_alt - count of tumour alt reads at a site
#' @param tumour_ref - count of tumour ref reads at a site
#' @param host_alt - count of host alt reads at a site
#' @param host_ref - count of host ref reads at a site
log_odds_ratio_variance <- function(tumour_alt, tumour_ref, host_alt, host_ref) {
    1 / (tumour_alt + 0.5) + 1 / (tumour_ref + 0.5) + 1 / (host_alt + 0.5) + 1 / (host_ref + 0.5)
}
