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

#' Reflect values in v around 0.5
#' (assumes 0 <= val <= 1, for all val in v)
#' @param v A vector of numbers in the interval [0, 1]
#' @examples
#' mirror(c(0.1, 0.2, 0.6, 0.7, 0.4, 0.8))
#' # [1] 0.1 0.2 0.4 0.3 0.4 0.2
#' @export
mirror <- function(v) ifelse(v > 0.5, 1 - v, v)

#' Clamp values in \code{v} to \code{lower} and \code{upper}, so that
#' any value \code{val} in \code{v} for which \code{val}
#' < \code{lower} is replaced by \code{lower}, and
#' \code{val} > \code{upper} is replaced by \code{upper}
#' @param v A vector of numbers.
#' @param lower A number. Any values in `v` below this are replaced.
#' @param upper A number. Any values in `v` above this are replaced.
#' @return A vector like \code{v}, with all values less than
#' \code{lower} or greater than \code{upper} replaced.
#' @examples
#' clamp(c(0, 0.4, 0.6, 1.0), 0.001, 0.999)
#' # [1] 0.001 0.400 0.600 0.999
#' @export
clamp <- function(v, lower, upper) {
    pmin(pmax(v, lower), upper)
}

scale_to_unit <- function(v) {
    r <- range(v)
    new <- (v - r[1]) / (r[2] - r[1])
    list(vals=new, undo_params=r)
}

undo_scale_to_unit <- function(v, undo_params) {
    v * (undo_params[2] - undo_params[1]) + undo_params[1]
}

#' Calculate the expected VAF for a SNP at the given purity, ploidy and logR,
#' if that SNP were LOH or otherwise homozygous in the founder
#' Can be used produce a segment threshold, above which a SNP is likely to have been
#' HET in the founder.
#' @param purity Estimated tumour purity (real value [0..1])
#' @param logr (numeric) Log (base-2) of the read depth ratio between tumour and normal
#' @param ploidy (numeric) Estimated overall tumour ploidy
#' @return The estimated VAF, had the original state of the SNP been homozygous
#' (real value [0..1]). NB: This is hypothetical, not a prediction.
#' @examples
#' expected_loh_vaf(0.817, c(1.713900, 1.607939, 1.629899, 1.570724, 1.666647), 2.114)
#' # [1] 0.02665127 0.02868239 0.02824911 0.02943189 0.02753864
#' @export
expected_loh_vaf <- function(purity, logr, ploidy) {
    #(1 - purity) / ((1 - purity) + 2^logr * (2 * (1 - purity) + purity * ploidy))
    totalcn <- calc_total_copynumber(logr, purity, ploidy)
    (1 - purity) / (pmax(0, totalcn) * purity + 2 * (1 - purity))
}

#' Returns segmented logR using pcf from the copynumber package
#' Intended for finishing off already segmented data, so ideally
#' the logr vector will be fairly short, and the penalty fairly
#' high
#' @importFrom "copynumber" pcf
#' @export
#' @param chr (character) Chromosome name
#' @param pos (integer vector) Genome positions
#' @param logr (real vector) logR values at each position
#' @param penalty (real) amount to penalise the PCF algorithm
#' creating a new segment (higher values -> fewer segments)
#' @param kmin (int) minimum number of data points in aa segment
resegment_logr <- function(chr, pos, logr, penalty=100, kmin=5) {
    df <- data.frame(CHROM = chr, POS = pos, LOGR = logr)
    suppressWarnings(segments <- pcf(df, gamma = penalty, kmin = kmin))
    setDT(segments)
    setnames(segments, old = c("chrom", "start.pos", "end.pos", "n.probes"), new = c("chr", "startpos", "endpos", "Nsnps"))
    segments[, width := endpos - startpos + 1]
    setcolorder(segments, c("chr", "startpos", "endpos", "width", "Nsnps"))
    segments[, .SD, .SDcols = c("chr", "startpos", "endpos", "width", "Nsnps")]
}

#' @export
listappend <- function(lst, item) {
    lst[[length(lst) + 1]] <- item
    return (lst)
}

