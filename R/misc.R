#' Rough estimate of purity
#' @export
estimate_purity <- function(host_vaf, tumour_vaf, logr,
                            hhetl=0.3, hhetu=0.7,
                            thoml=0.25, thomu=0.75,
                            lrl=-0.75, lru=0.75) {
    # approx_ploidy <- 2 / 2^mean(logr)

    mirror <- function(v) ifelse(v > 0.5, 1 - v, v)
    het_host_ix <- which(host_vaf > hhetl & host_vaf < hhetu)
    hom_tum_ix  <- which(tumour_vaf < thoml | tumour_vaf > thomu)
    tum_cn2_ix  <- which(logr > lrl & logr < lru)
    ix <- intersect(het_host_ix, intersect(hom_tum_ix, tum_cn2_ix))

    approx_purity <- 1 - 2 * mean(mirror(tumour_vaf[ix]))
    approx_ploidy <- 2 / (2^(mean(logr[tum_cn2_ix])) * approx_purity) - 2 / approx_purity + 2
    list(purity = approx_purity, ploidy = approx_ploidy)
}
