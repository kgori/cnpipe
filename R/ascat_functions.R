#' @importFrom "data.table" fread setcolorder
#' @export
process_tsv <- function(filename, min_host = 10, min_tumour = 1) {
    # Load data and select columns
    dt <- fread(filename, showProgress = FALSE)
    cols.to.keep <- c(grep("CHROM|POS|REF|ALT", colnames(dt)), 7, 8, 11, 12, 13, 14)
    baflogr <- dt[, cols.to.keep, with = FALSE]
    colnames(baflogr) <- c("chr", "cumpos", "REF", "ALT", "h_total", "h_alt", "t_total", "t_alt", "baf", "logr")
    baflogr[, "h_ref" := h_total - h_alt]
    baflogr[, "t_ref" := t_total - t_alt]
    setcolorder(baflogr, c("chr", "cumpos", "REF", "ALT", "h_ref", "h_alt", "t_ref", "t_alt", "baf", "logr", "h_total", "t_total"))

    # Filter low coverage SNVs
    mask <- baflogr[, h_total >= min_host & t_total >= min_tumour]

    # Remove unnecessary columns
    baflogr[, h_total := NULL]
    baflogr[, t_total := NULL]

    # Return filtered
    baflogr[mask]
}

calc_logr <- function(h_ref, h_alt, t_ref, t_alt) {
    h_tot <- h_ref + h_alt
    t_tot <- t_ref + t_alt
    log2( (t_tot / h_tot) / median((t_tot / h_tot), na.rm = T) )
}

#' @export
calc_total_copynumber <- function(logr, purity, ploidy) {
    ((2 * (1 - purity) + purity * ploidy) * 2^logr - 2 * (1 - purity)) / purity
}

#' Calculates host genotype
#' host = hom if (coverage ≥ 10, as above) and 0 reads ref | alt
#' host = het if ≥ 3 reads for both ALT and REF & .25 ≤ BAF ≤ .75
#' @param data = data frame / table with columns 'h_ref', 'h_alt', 't_ref'
#' and 't_alt'
#' @export
dep_calc_combined_genotype <- function(h_ref, h_alt, t_ref, t_alt,
                                   het_baf_interval = c(0.25, 0.75),
                                   hom_baf_boundaries = c(0.1, 0.9),
                                   min_reads = 3) {
    h_baf <- h_alt / (h_ref + h_alt)

    combined_gt <- ifelse(h_ref < min_reads & h_baf > hom_baf_boundaries[2],
                          ifelse(t_ref >= min_reads, "BB/A*", "BB/B"),
                          ifelse(h_alt < min_reads & h_baf < hom_baf_boundaries[1],
                                 ifelse(t_alt >= min_reads, "AA/B*", "AA/A"),
                                 ifelse(h_ref >= min_reads & h_alt >= min_reads & h_baf %between% het_baf_interval,
                                        "AB/*",
                                        NA)))
    combined_gt
}

#' @export
calc_host_genotype <- function(h_baf, groups = 3, min_reads = 3, sample_size = 10000, tolerance = 1e-03, plot = FALSE, estimate_uncertainty = FALSE) {
    mixture_data <- h_baf[!is.na(h_baf) & h_baf > 0 & h_baf < 1]
    mixture_data_sample <- sample(mixture_data, sample_size, replace = sample_size > length(mixture_data))
    inits <- get_init(mixture_data, groups)
    # print(inits)
    mmod <- estimate_mixture(mixture_data_sample, inits$params, inits$pi, tol = tolerance, plot = plot)

    # Don't assume means are in sorted order (but they should be)
    means <- mmod$means
    homref_class <- which.min(means)
    homalt_class <- which.max(means)
    het_class <- setdiff(1:3, c(homref_class, homalt_class))

    # Apply prediction to full data sample
    weights <- get_weights(h_baf, mmod$params, mmod$pi)$weights
    classification <- apply(weights, 1, which.max)

    genotype <- ifelse(classification == homref_class,
                       "AA",
                       ifelse(classification == het_class,
                              "AB",
                              ifelse(classification == homalt_class,
                                     "BB",
                                     NA)))

    if (estimate_uncertainty) {
        uncertainty <- 1 - apply(weights, 1, max)
        return (list(genotype = genotype, uncertainty = uncertainty))
    } else {
        return (genotype)
    }
}

#' @export
calc_host_genotype_quick <- function(h_baf, groups = 3, plot = FALSE) {
    dens <- density(c(-h_baf, h_baf, 1+h_baf), from = 0, to = 1)
    homref_bound <- dens$x[dens$x < 0.5][which.min(dens$y[dens$x < 0.5])]
    homalt_bound <- dens$x[dens$x > 0.5][which.min(dens$y[dens$x > 0.5])]
    print(c(homref_bound, homalt_bound))
    if (plot) {
        hist(h_baf, breaks = seq(0, 1, length.out = 41),
             border = "black", col = "paleturquoise1",
             probability = TRUE, xlim = c(0, 1))#, ylim = c(0, max(dens$y)))
        lines(dens, lwd = 2, col = "steelblue4", xlim = c(0, 1))
        abline(v = c(homalt_bound, homref_bound), lty = 2, lwd = 2, col = "palevioletred")
    }
    ifelse(h_baf < homref_bound, "AA",
           ifelse(h_baf < homalt_bound, "AB",
                  "BB"))
}

#' @export
calc_combined_genotype <- function(host_genotype, t_ref, t_alt, min_reads = 3) {
    ifelse(host_genotype == "AB",
           "AB/*",
           ifelse(host_genotype == "AA",
                  ifelse(t_alt >= min_reads, "AA/B*", "AA/A"),
                  ifelse(host_genotype == "BB",
                         ifelse(t_ref >= min_reads, "BB/A*", "BB/B"),
                         NA)))
}

#' @export
calc_allele_specific_copynumber <- function(htgeno, baf, logr, purity, ploidy) {
    na <- ifelse(htgeno == "AB/*",
                 # Host is heterozygous:
                 tumour_a_allele_copy_number(logr = logr, baf = baf,
                                             purity = purity, tumour_ploidy = ploidy,
                                             ha = 1, hb = 1),

                 # else:
                 ifelse(htgeno %in% c("AA/B*", "AA/A"),
                        # Host is homozygous ref:
                        tumour_a_allele_copy_number(logr = logr, baf = baf,
                                                    purity = purity, tumour_ploidy = ploidy,
                                                    ha = 2, hb = 0),
                        # else:
                        ifelse(htgeno %in% c("BB/A*", "BB/B"),
                               # Host is homozygous alt:
                               tumour_a_allele_copy_number(logr = logr, baf = baf,
                                                           purity = purity, tumour_ploidy = ploidy,
                                                           ha = 0, hb = 2),

                               # Else ???:
                               NA)))

    nb <- ifelse(htgeno == "AB/*",
                 # Host is heterozygous:
                 tumour_b_allele_copy_number(logr = logr, baf = baf,
                                             purity = purity, tumour_ploidy = ploidy,
                                             ha = 1, hb = 1),
                 ifelse(htgeno %in% c("AA/B*", "AA/A"),
                        # Host is homozygous ref:
                        tumour_b_allele_copy_number(logr = logr, baf = baf,
                                                    purity = purity, tumour_ploidy = ploidy,
                                                    ha = 2, hb = 0),

                        ifelse(htgeno %in% c("BB/A*", "BB/B"),
                               # Host is homozygous alt:
                               tumour_b_allele_copy_number(logr = logr, baf = baf,
                                                           purity = purity, tumour_ploidy = ploidy,
                                                           ha = 0, hb = 2),

                               # Else ???:
                               NA)))

    list(na = na, nb = nb)
}

#' Convert vectors of start and end points of segments to
#' vector of indices of points within those segments
#' Inverse of index_to_segments
#' @export
segments_to_index <- function(starts, ends) {
    unlist(mapply(starts, ends, FUN = ':'))
}

#' Convert vector of indices of points within segments to a list of
#' vectors of start and end points of those segments
#' Inverse of segments_to_index
#' @export
index_to_segments <- function(v) {
    r <- rle(v - 1:length(v))
    cs <- cumsum(r$lengths)
    starts <- v[c(1, cs[1:(length(cs)-1)] + 1)]
    ends <- v[cs]
    list(starts=starts, ends=ends)
}

#' @importFrom "stats" density
#' @export
findpeaks <- function (x, thresh = 0.00001, from = -.75, to = 1.75, mindens = .25) {
    xdens <- density(x, kernel = "gaussian", from = from, to = to)
    pks <- which(diff(sign(diff(xdens$y, na.pad = FALSE)), na.pad = FALSE) < 0) + 2
    pks <- pks[xdens$y[pks - 1] - xdens$y[pks] > thresh & xdens$y[pks] > mindens]
    return(xdens$x[pks])
}

# #' Segment copy number using an HMM. Separates a population of points with a "small" value
# #' from everything with a "larger" value, e.g. all copynumber==0 SNVs from all copynumber>=1 SNVs
# #' @param cn Copy number data
# #' @param prob01 Probability of transitioning from state CN0 to state CN1+
# #' @param prob10 Probability of transitioning from state CN1+ to state CN0
# #' @param mean0 Mean of state 0
# #' @param mean1 Mean of state 1
# #' @param python_exe Full path to python executable
# #' @param params String describing the parameters hmmlearn should optimise -
# #' s = start probs, m = means, c = covars, t = transition probs
# #' @param return_smaller_state Invert the selection, i.e. return the set of "large" values
# #' @param plot Histogram of data points belonging to state 0 (blue) or state 1 (red)
# #' @return list of segment starts, segment ends, and the trained model
# #' @importFrom "reticulate" use_python source_python
# #' @export
# hmm_segmentation <- function(cn, prob01, prob10, mean0, mean1, minrun = 1000, python_exe,
#                              params = "smc", return_smaller_state = TRUE, plot = FALSE) {
#     use_python(python_exe)
#     modelfile <- system.file("python", "gaussianhmm.py", package = "cnpipe")
#     stopifnot(file.exists(modelfile))
#     source_python(modelfile)
#
#     model <- model_builder(prob_change_state_01 = prob01, prob_change_state_10 = prob10,
#                            mean_state0 = mean0, mean_state1 = mean1, params = params)
#     data <- matrix(cn, ncol = 1)
#     model$fit(data)
#     prediction <- as.vector(model$predict(data))
#     runlengths <- rle(prediction)
#     smaller_state <- ifelse(model$means_[1] < model$means_[2], 0, 1)
#     larger_state <- ifelse(model$means_[1] < model$means_[2], 1, 0)
#     if (return_smaller_state) {
#         runlengths$values[runlengths$lengths < minrun] <- larger_state
#     }  else {
#         runlengths$values[runlengths$lengths < minrun] <- smaller_state
#     }
#
#     if (plot) {
#         p1 <- hist(cn[prediction==0], breaks = 100)
#         p2 <- hist(cn[prediction==1], breaks = 100)
#         lim <- range(cn)
#         plot(p1, xlim = lim, col = rgb(0,0,1,1/4), border = NA, main = NA, xlab = "NMIN", freq = FALSE)
#         plot(p2, xlim = lim, col = rgb(1,0,0,1/4), border = NA, add = T, freq = FALSE)
#
#     }
#
#     if(return_smaller_state) {
#         if (smaller_state == 0) {
#             prediction <- !inverse.rle(runlengths)
#         } else {
#             prediction <- !!inverse.rle(runlengths)
#         }
#     } else {
#         if (larger_state == 0) {
#             prediction <- !inverse.rle(runlengths)
#         } else {
#             prediction <- !!inverse.rle(runlengths)
#         }
#     }
#     starts <- which(c(prediction, F) & !c(F, prediction))
#     ends <- which(!c(prediction, F) & c(F, prediction))
#     list(starts = starts, ends = ends, model = model)
# }

#' Write the 4 required files for ascat
#' @param data - contains data needed by ascat
#' @param filenames - list of 4 filenames, accessed as: tumour_logr, normal_logr, tumour_baf, normal_baf
#' @export
write_ascat_data <- function(data, filenames) {
    stopifnot(length(filenames) == 4)
    fwrite(x = data.frame(chr = data$chr,
                          pos = data$cumpos,
                          'tumour' = data$logr),
           file = filenames$tumour_logr,
           sep = "\t",
           row.names = TRUE)

    fwrite(x = data.frame(chr = data$chr,
                          pos = data$cumpos,
                          'normal' = 0),
           file = filenames$normal_logr,
           sep = "\t",
           row.names = TRUE)

    fwrite(x = data.frame(chr = data$chr,
                          pos = data$cumpos,
                          'tumour' = abs(data$baf - sample(x = c(0,1), size = nrow(data), replace = T))),
           file = filenames$tumour_baf,
           sep = "\t",
           row.names = TRUE)

    fwrite(x = data.frame(chr = data$chr,
                          pos = data$cumpos,
                          'normal' = abs(data$h_alt/(data$h_ref+data$h_alt) - sample(x = c(0,1), size = nrow(data), replace = TRUE))),
           file = filenames$normal_baf,
           sep = "\t",
           row.names = TRUE)
}

#' @export
tumour_a_allele_copy_number <- function(logr, baf, purity, tumour_ploidy, ha = 1, hb = 1) {
    (2^logr * (1 - baf) * ((ha + hb) * (1 - purity) + purity * tumour_ploidy) - ha * (1 - purity)) / purity
}

#' @export
tumour_b_allele_copy_number <- function(logr, baf, purity, tumour_ploidy, ha = 1, hb = 1) {
    (2^logr * baf * ((ha + hb) * (1 - purity) + purity * tumour_ploidy) - hb * (1 - purity)) / purity
}

#' @export
tumour_total_copy_number <- function(logr, purity, tumour_ploidy, ha = 1, hb = 1) {
    (2^logr * ((ha + hb) * (1 - purity) + purity * tumour_ploidy) - (ha + hb) * (1 - purity)) / purity
}

#' Convert an ascat object (not an ascat result) to a data table
#' containing SNP positions, logR and BAF at each position,
#' and segmented average LogR and BAF where available
#' @importFrom "data.table" as.data.table copy foverlaps setnames setkey setorder setcolorder
#' @param ascat_obj ASCAT data object, including segmentation information
#' @param ascat_result ASCAT analysis result object
#' @return data table representation of ascat data
#' @export
ascatobj_to_datatable <- function(ascat_obj, ascat_result = NULL) {
    ascat_obj <- copy(ascat_obj)
    ascat_result <- copy(ascat_result)

    # Tidying up object names
    colnames(ascat_obj$Tumor_BAF) <- sub("\\..*$", "", colnames(ascat_obj$Tumor_BAF))
    colnames(ascat_obj$Tumor_LogR) <- sub("\\..*$", "", colnames(ascat_obj$Tumor_LogR))
    colnames(ascat_obj$Tumor_LogR_segmented) <- sub("\\..*$", "", colnames(ascat_obj$Tumor_LogR_segmented))
    ascat_obj$samples <- sub("\\..*$", "", ascat_obj$samples)

    snp_pos <- as.data.table(ascat_obj$SNPpos)
    snp_pos[, ix := .I]
    i <- 0
    dts <- list()

    for (samplename in ascat_obj$samples) {
        i <- i + 1

        base_dt <- as.data.table(ascat_obj$SNPpos)
        base_dt[, chrom := as.character(chrom)]
        base_dt[, snp_index := .I]
        base_dt[, sample := samplename]
        base_dt[, ("logR") := ascat_obj$Tumor_LogR[, samplename]]
        base_dt[, ("segmented_logR") := ascat_obj$Tumor_LogR_segmented[, samplename]]
        base_dt[, ("vaf") := ascat_obj$Tumor_BAF[, samplename]]
        segbaf <- ascat_obj$Tumor_BAF_segmented[[i]]
        rn <- rownames(segbaf)
        base_dt[snp_index %in% rn, ("segmented_vaf") := segbaf[, 1]]
        dts[[samplename]] <- base_dt
    }

    dts <- rbindlist(dts)

    if (!is.null(ascat_result)) {
        dts[, c("start.pos", "end.pos") := POS]
        setkey(dts, sample, chrom, start.pos, end.pos)

        seg <- copy(ascat_result$segments_raw)
        setDT(seg)
        seg[, sampleSegID := .I]
        seg[, sample := sub("\\..*$", "", sample)]
        setnames(seg, old = c("chr", "startpos", "endpos"), new = c("chrom", "start.pos", "end.pos"))
        setkey(seg, sample, chrom, start.pos, end.pos)

        dts <- foverlaps(dts, seg)
        dts[, datasetSegID := paste0(unique(sampleSegID), collapse="_"), .(chrom, POS)]
        dts[, datasetSegID := as.numeric(as.factor(datasetSegID))]

        dts[, c("start.pos", "end.pos", "i.start.pos", "i.end.pos") := NULL]
        dts[, c("orig.start.pos", "orig.end.pos", "orig.nSNPs") := .(min(POS), max(POS), .N), by = .(sample, sampleSegID)]
        dts[, c("start.pos", "end.pos", "nSNPs") := .(min(POS), max(POS), .N), by = .(sample, datasetSegID)]
        setcolorder(dts, c("snp_index", "chrom", "POS", "sample", "nMajor", "nMinor","nAraw", "nBraw", "logR", "vaf", "segmented_logR", "segmented_vaf",
                           "sampleSegID", "datasetSegID", "start.pos", "end.pos", "nSNPs", "orig.start.pos", "orig.end.pos", "orig.nSNPs"))
        setorder(dts, snp_index, chrom, POS, sample)
    }
    return (dts)
}

#' @export
idealised_vaf <- function(ta, tb, ha, hb, purity) {
    (hb * (1 - purity) + tb * purity) / ((ha+hb) * (1-purity) + (ta+tb) * purity)
}

#' @export
idealised_logr <- function(ta, tb, ha, hb, purity, ploidy, host_ploidy = 2) {
    log2(((ha+hb) * (1-purity) + (ta+tb) * purity) / (((1-purity) * host_ploidy + purity * ploidy)))
}

replace_na <- function(val, default) {
    ifelse(is.na(val), default, val)
}

#' Simulate VAF according to certain parameters
#' @export
simulated_vaf <- function(samplesize, depth, ta, tb, ha, hb, purity, tploidy = 2, hploidy = 2, balance = 0) {
    stopifnot(purity >= 0 & purity <= 1)

    # Adjust the sequencing depth of this region from the average `depth`
    # based on this region's total copy number and the tumour and host average ploidies
    region_depth <- (((ta+tb) / tploidy) * purity + (1 - purity) * ((ha+hb) / hploidy)) * depth

    # Work out the expected number of reads based on this region's depth
    expected_ta <- replace_na(region_depth * purity * ta / (ta+tb), 0)
    expected_tb <- replace_na(region_depth * purity * tb / (ta+tb), 0)
    expected_ha <- replace_na(region_depth * (1-purity) * ha / (ha+hb), 0)
    expected_hb <- replace_na(region_depth * (1-purity) * hb / (ha+hb), 0)

    # Sample read counts from a Poisson distribution
    nta <- rpois(samplesize, expected_ta + 0.01)
    ntb <- rpois(samplesize, expected_tb + 0.01)
    nha <- rpois(samplesize, expected_ha + 0.01)
    nhb <- rpois(samplesize, expected_hb + 0.01)

    vaf <- ((ntb+nhb) / (nta+ntb+nha+nhb))
    expected_vaf <- ((expected_tb+expected_hb) / (expected_ta+expected_tb+expected_ha+expected_hb))

    # Rebalance the SNPs onto major and minor alleles
    ix <- seq(1, balance * length(vaf))
    vaf[ix] <- 1 - vaf[ix]
    list(data = vaf, expectation = expected_vaf)
}

#' Simulate logr according to certain parameters
#' @export
simulated_logr <- function(samplesize, depth, tn, hn, purity, tploidy = 2, hploidy = 2) {
    stopifnot(purity >= 0 & purity <= 1)

    # Work out the expected number of reads in this segment in a host sample
    # sequenced to equal depth
    expected_h <- (hn / hploidy) * depth

    # Work out the expected number of reads mapping to this segment in the tumour sample,
    # given the average depth, ploidy and purity
    expected_t <- ((tn / tploidy) * purity + (1 - purity) * (hn / hploidy)) * depth

    # Sample read counts from a Poisson distribution
    nt <- rpois(samplesize, expected_t + 0.01)
    nh <- rpois(samplesize, expected_h + 0.01)
    logr <- log2(nt/nh)
    expected_logr <- log2(expected_t / expected_h)
    list(data = logr, expectation = expected_logr)
}


#' Annotate `breaks` (=ascat_result$segments) so each segment is given two IDs:
#'   dataset_segment_id - uniquely identifies the segment (chr, startpos, endpos) within the entire dataset
#'                      - (a segment can occur in multiple samples)
#'   sample_segment_id - gives the segment an index within the sample
#' And two frequencies:
#'   breakpoint_frequency - number of times a given startpos is observed
#'   segment_frequency - number of times a unique segment is observed
#' @export
annotate_ascat_segments <- function(breaks) {

    # Assign sample segment ids
    setorder(breaks, sample, chr, startpos, endpos)
    for (samplename in unique(breaks$sample)) {
        breaks[sample==samplename, sample_segment_id := seq(.N)]
    }

    # Count dataset segment and breakpoint frequencies
    breaks[, breakpoint_frequency := .N, by = .(chr, startpos)]
    breaks[, segment_frequency := .N, by = .(chr, startpos, endpos)]

    # Collect the unique segments (can occur in multiple samples)
    ubk <- unique(breaks, by = c("chr", "startpos", "endpos"))
    setorder(ubk, sample, chr, startpos, endpos)
    ubk[, dataset_segment_id := .I]

    # Transfer the unique IDs to the breakpoints list
    setkey(ubk, chr, startpos, endpos)
    setkey(breaks, chr, startpos, endpos)
    breaks <- breaks[ubk][, -c("i.sample", "i.sample_segment_id", "i.nMajor", "i.nMinor",
                               "i.breakpoint_frequency", "i.segment_frequency")][order(sample, chr, startpos, endpos)]
    setcolorder(breaks, c("sample", "chr", "startpos", "endpos", "nMajor", "nMinor", "sample_segment_id", "dataset_segment_id", "segment_frequency"))
    setkey(breaks, sample, chr, startpos, endpos)
    return (breaks)
}


#' Create a new table of segments that inlcudes every breakpoint from every sample in
#' `datatable`, where datatable = ascat_result$segments
#' @importFrom "stringr" str_sort
#' @export
disjoin_ascat_segments <- function(datatable) {
    .disjoin.chr <- function(dt) {
        # data.table disjoin
        starts <- unique(dt$startpos)
        ends <- unique(dt$endpos)
        adjstart <- head(sort(unique(c(starts, ends + 1))), -1)
        adjend <- tail(sort(unique(c(ends, starts - 1))), -1)
        adj <- data.table(startpos=adjstart, endpos=adjend, width=adjend - adjstart + 1)
        setkey(adj, startpos, endpos)
        unique(foverlaps(dt, adj, nomatch=0L, minoverlap=1L),
               by = c("startpos", "endpos"))[, -c("i.startpos", "i.endpos")]
    }
    res <- datatable[, .disjoin.chr(.SD), by = chr, .SDcols = c("startpos", "endpos")]

    # Set sort order
    chroms <- stringr::str_sort(unique(res$chr))
    res[, chr := factor(chr, levels = chroms)]
    setorder(res, chr, startpos, endpos)
    res[, chr := as.character(chr)]
    res$ID <- res[, .I]
    setkey(res, chr, startpos, endpos)
    return(res)
}

#' Post-process a disjointed segments list to remove segments that cover 0 SNPs
#' Adds a column, Nsnps, the number of SNPs covered by each segment
#' @param snp_list (data.table) A list of SNPs with columns chr (chromosome),
#' startpos (SNP position) and endpos (=startpos). Other columns can be present,
#' but ignored.
#' @param segments (data.table) A list of segments with columns chr (chromosome),
#' startpos and endpos (start and end positions of the segments)
#' @importFrom "data.table" setkey key foverlaps
#' @export
filter_segments_by_snplist <- function(snp_list, segments) {
    has_cols <- function(dt, names) {
        all(intersect(colnames(dt), names) == names)
    }

    stopifnot(has_cols(snp_list, c("chr", "startpos", "endpos")))
    stopifnot(has_cols(segments, c("chr", "startpos", "endpos")))

    if(!has_cols(segments, c("width"))) {
        segments[, width := endpos - startpos + 1]
    }

    setkey(segments, chr, startpos, endpos)
    do.call(setkey, append(list(snp_list), as.list(key(segments))))

    foverlaps(snp_list, segments, nomatch=0L)[, .(Nsnps = .N), by = .(chr, startpos, endpos, width)]
}
