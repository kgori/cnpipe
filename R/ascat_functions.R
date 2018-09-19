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
calc_host_genotype <- function(h_ref, h_alt, groups = 3, min_reads = 3, sample_size = 10000, plot = FALSE) {
    h_baf <- h_alt / (h_ref + h_alt)
    mixture_data <- h_baf[h_baf > 0 & h_baf < 1]
    mixture_data_sample <- sample(mixture_data, sample_size, replace = sample_size > length(mixture_data))
    inits <- get_init(mixture_data, groups)
    mmod <- estimate_mixture(mixture_data_sample, inits$params, inits$pi, tol = 1e-03, plot = plot)

    if (plot) {
        lines(mmod$data[order(mmod$data)], mmod$uncertainty[order(mmod$data)])
        lower <- beta_boundary(mmod$params[1, 1], mmod$params[1, 2], mmod$params[2, 1], mmod$params[2, 2], mmod$pi[1], mmod$pi[2], log = TRUE)
        if (groups == 3) {
            upper <- beta_boundary(mmod$params[2, 1], mmod$params[2, 2], mmod$params[3, 1], mmod$params[3, 2], mmod$pi[2], mmod$pi[3], log = TRUE)
            abline(v = c(lower, upper), col = "red", lty = 2)
        } else {
            abline(v = lower, col = "red", lty = 2)
        }
    }

    # Don't assume means are in sorted order (but they should be)
    means <- mmod$means
    homref_class <- which.min(means)
    homalt_class <- which.max(means)
    het_class <- setdiff(1:3, c(homref_class, homalt_class))

    # Apply prediction to full data sample
    weights <- get_weights(h_baf, mmod$params, mmod$pi)$weights
    classification <- apply(weights, 1, which.max)

    ifelse(classification == homref_class,
           "AA",
           ifelse(classification == het_class,
                  "AB",
                  ifelse(classification == homalt_class,
                         "BB",
                         NA)))
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

# @importFrom "matrixStats" logSumExp
gaussian_boundary <- function(mean1, mean2, sd1, sd2, log=TRUE) {
    x <- seq(mean1, mean2, length.out = 1001)
    df <- data.frame(dens1 = dnorm(x, mean1, sd1, log=log),
                     dens2 = dnorm(x, mean2, sd2, log=log))
    y <- apply(df, 1, function(row) abs(row[1] - row[2]))
    result <- x[which.min(y)]
    # plot(x, y, type = "l", ylim = c(min(c(y, df$dens1, df$dens2)), max(y)))
    # lines(x, df$dens1, lty=2)
    # lines(x, df$dens2, lty=2)
    abline(v = result, col = "red")
    result
}

#' @export
calc_allele_specific_copynumber <- function(htgeno, baf, logr, purity, ploidy) {
    na <- ifelse(htgeno == "AB/*",
                 # Host is heterozygous:
                 ((2 * (1 - purity) + purity * ploidy) * (1 - baf) * 2^logr - (1 - purity)) / purity,

                 # else:
                 ifelse(htgeno %in% c("AA/B*", "AA/A"),
                        # Host is homozygous ref:
                        ((2 * (1 - purity) + purity * ploidy) * (1 - baf) * 2^logr - 2 * (1 - purity)) / purity,
                        # else:
                        ifelse(htgeno %in% c("BB/A*", "BB/B"),
                               # Host is homozygous alt:
                               ((2 * (1 - purity) + purity * ploidy) * (1 - baf) * 2^logr) / purity,

                               # Else ???:
                               NA)))

    nb <- ifelse(htgeno == "AB/*",
                 # Host is heterozygous:
                 ((2 * (1 - purity) + purity * ploidy) * baf * 2^logr - (1 - purity)) / purity,
                 ifelse(htgeno %in% c("AA/B*", "AA/A*"),
                        # Host is homozygous ref:
                        ((2 * (1 - purity) + purity * ploidy) * baf * 2^logr) / purity,

                        ifelse(htgeno %in% c("BB/A*", "BB/B*"),
                               # Host is homozygous alt:
                               ((2 * (1 - purity) + purity * ploidy) * baf * 2^logr - 2 * (1 - purity)) / purity,

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

#' Segment copy number using an HMM. Separates a population of points with a "small" value
#' from everything with a "larger" value, e.g. all copynumber==0 SNVs from all copynumber>=1 SNVs
#' @param cn Copy number data
#' @param prob01 Probability of transitioning from state CN0 to state CN1+
#' @param prob10 Probability of transitioning from state CN1+ to state CN0
#' @param mean0 Mean of state 0
#' @param mean1 Mean of state 1
#' @param python_exe Full path to python executable
#' @param params String describing the parameters hmmlearn should optimise -
#' s = start probs, m = means, c = covars, t = transition probs
#' @param return_smaller_state Invert the selection, i.e. return the set of "large" values
#' @param plot Histogram of data points belonging to state 0 (blue) or state 1 (red)
#' @return list of segment starts, segment ends, and the trained model
#' @importFrom "reticulate" use_python source_python
#' @export
hmm_segmentation <- function(cn, prob01, prob10, mean0, mean1, minrun = 1000, python_exe,
                             params = "smc", return_smaller_state = TRUE, plot = FALSE) {
    use_python(python_exe)
    modelfile <- system.file("python", "gaussianhmm.py", package = "cnpipe")
    stopifnot(file.exists(modelfile))
    source_python(modelfile)

    model <- model_builder(prob_change_state_01 = prob01, prob_change_state_10 = prob10,
                           mean_state0 = mean0, mean_state1 = mean1, params = params)
    data <- matrix(cn, ncol = 1)
    model$fit(data)
    prediction <- as.vector(model$predict(data))
    runlengths <- rle(prediction)
    smaller_state <- ifelse(model$means_[1] < model$means_[2], 0, 1)
    larger_state <- ifelse(model$means_[1] < model$means_[2], 1, 0)
    if (return_smaller_state) {
        runlengths$values[runlengths$lengths < minrun] <- larger_state
    }  else {
        runlengths$values[runlengths$lengths < minrun] <- smaller_state
    }

    if (plot) {
        p1 <- hist(cn[prediction==0], breaks = 100)
        p2 <- hist(cn[prediction==1], breaks = 100)
        lim <- range(cn)
        plot(p1, xlim = lim, col = rgb(0,0,1,1/4), border = NA, main = NA, xlab = "NMIN", freq = FALSE)
        plot(p2, xlim = lim, col = rgb(1,0,0,1/4), border = NA, add = T, freq = FALSE)

    }

    if(return_smaller_state) {
        if (smaller_state == 0) {
            prediction <- !inverse.rle(runlengths)
        } else {
            prediction <- !!inverse.rle(runlengths)
        }
    } else {
        if (larger_state == 0) {
            prediction <- !inverse.rle(runlengths)
        } else {
            prediction <- !!inverse.rle(runlengths)
        }
    }
    starts <- which(c(prediction, F) & !c(F, prediction))
    ends <- which(!c(prediction, F) & c(F, prediction))
    list(starts = starts, ends = ends, model = model)
}

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
