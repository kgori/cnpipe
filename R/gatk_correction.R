#' Clip numeric vector so that its lowest element is >= low,
#' and its highest is <= high
clip <- function(v, low, high) {
    pmax(low, pmin(v, high))
}

safelog2 <- function(x) {
    clip(log2(x), -31, 31)
}

#' Creates an object to do GATK logR standardization and denoising using a panel of normals.
#' This is simplified compared to the GATK procedure, as it doesn't apply any filters to the
#' data, except for an optional winsorization. As such, it is advisable to filter out rows of
#' the host matrix that have very low or high coverage, or which have row-median of zero.
#' @importFrom "R6" R6Class
#' @importFrom "matrixStats" rowMedians
#' @export
logRCorrector = R6Class(
    "logRCorrector",
    public=list(
        #' @field row_medians Row medians of the input host read count matrix
        row_medians=NA,

        #' @field U Result of Singular Value Decomposition of processed host read counts,
        #' used for denoising tumour logR.
        U=NA,

        #' @description
        #' Constructor. Builds new object using host read count matrix as input
        #' @param host_matrix Matrix of read counts, with samples as columns.
        #' @param winsorization_quantile Values outside the quantile range set here
        #' (q, 1-q) are clamped to the values at those quantiles.
        #' Can be set to 0 to disable.
        initialize = function(host_matrix, winsorization_quantile=0.001) {
            stopifnot(winsorization_quantile >= 0 & winsorization_quantile <= 1)

            host_matrix <-  sweep(host_matrix, 2, colMedians(host_matrix), "/")
            rm <- rowMedians(host_matrix)
            rm[rm==0]<-1
            host_matrix <- sweep(host_matrix, 1, rm, "/")
            host_matrix <- apply(host_matrix, 2, function(col) clip(col,
                                                                    quantile(col, winsorization_quantile),
                                                                    quantile(col, 1-winsorization_quantile)))
            ix=which(host_matrix==0, arr.ind=TRUE)
            host_matrix[ix] <- rowMedians(host_matrix)[ix[, "row"]]
            host_matrix <- sweep(host_matrix, 2, colMedians(host_matrix), "/")
            host_matrix <- apply(host_matrix, 2, safelog2)
            u <- svd(host_matrix)$u
            self$row_medians = rm
            self$U = u
        },

        #' @description
        #' Calculate standardized (not denoised) logR from tumour read counts
        #' @param tumour_matrix Matrix of tumour read counts, samples in columns
        standardize = function(tumour_matrix) {
            stopifnot(nrow(tumour_matrix) == nrow(self$U))
            tumour_matrix <- sweep(tumour_matrix, 2, colMedians(tumour_matrix), "/")
            tumour_matrix <- sweep(tumour_matrix, 1, self$row_medians, "/")
            tumour_matrix <- apply(tumour_matrix, 2, safelog2)
            tumour_matrix
        },

        #' @description
        #' Calculate denoised logR from tumour read counts.
        #' @param tumour_matrix Matrix of tumour read counts, samples in columns
        #' @param num_svs Number of singular values to use for denoising.
        #' @param randomized Set to TRUE to use a random selection of
        #' singular values. Set to FALSE to use the largest singular values.
        denoise = function(tumour_matrix, num_svs, randomized=FALSE) {
            stopifnot(nrow(tumour_matrix) == nrow(self$U))
            if (num_svs > ncol(self$U)) num_svs <- ncol(self$U)
            if (num_svs < 1) num_svs <- 1

            tumour_matrix <- self$standardize(tumour_matrix)
            if (randomized) {
                u_select <- self$U[, sample.int(ncol(.self$U), num_svs)]
            } else {
                u_select <- self$U[, 1:num_svs]
            }
            tumour_matrix - u_select %*% (t(u_select) %*% tumour_matrix)
        })
)


#' Obtain read count normalisation info from a panel of normals
#' @param data (data.table / data.frame) Table of read counts, with columns representing samples
#' and rows representing genomic locations or intervals
#' @param filter1 (real) Filter out rows where the median value is below the quantile threshold
#' given by this parameter (default = 0.1)
#' @param filter2 (real) Filter out rows where more than this threshold proportion of samples have
#' zero read coverage (default = 0.05)
#' @param filter3 (real) Clamp (Winsorize) values to the quantiles given by this value (and 1 - this value)
#' @param nsvs (int) Number of Singular Values to retain (cannot exceed number of samples)
#' @importFrom "matrixStats" rowMedians
#' @importFrom "logging" loginfo
#' @importFrom "data.table" setDT
#' @export
gatk_decompose <- function(data, filter1 = 0.1, filter2 = 0.05, filter3 = 0.001, nsvs = 20, verbose = TRUE) {
    setDT(data) # Must be a data.table, not a data.frame
    samples <- colnames(data)
    if (verbose) loginfo(paste0("Data dimensions:", nrow(data), " regions; ", ncol(data), " samples.\n"), logger = "LOGR")
    if (verbose) loginfo("Making sure data is numeric...", logger = "LOGR")
    data[, (samples) := lapply(.SD, function(value) as.numeric(value)), .SDcols = samples]

    # Secret undocumented GATK step 1A - convert read counts to fractional read counts
    if (verbose) loginfo("Converting to fractional read counts...", logger = "LOGR")
    data[, (samples) := lapply(.SD, function(value) value / sum(value)), .SDcols = samples]

    # GATK step 2 - Calculate target medians
    if (verbose) loginfo("Calculating median across all samples...", logger = "LOGR")
    rowmedian = data[, rowMedians(as.matrix(.SD)), .SDcols = samples]

    # GATK step 3 - Apply 1st filter - remove intervals with median across all samples below 10%ile
    if (verbose) loginfo(paste0("Filtering out regions with a median across all samples below the ", filter1*100, "%% quantile..."), logger = "LOGR")
    mask <- (rowmedian >= quantile(rowmedian, filter1) & rowmedian != 0) # filters out bad rows - a vector of length == nrow(data),
    # TRUE for each row to retain, FALSE for each row to remove

    # GATK step 4: Divide each row by row median
    rowmedian[rowmedian == 0] <- 1  # dummy count to avoid dividing by zero
    if (verbose) loginfo("Dividing each region by its median across all samples...", logger = "LOGR")
    data[mask] <- data[mask] / rowmedian[mask]

    # GATK step 5 - filter out samples with >5% 0-coverage targets (SKIP)

    # GATK step 6 - filter out targets with >5% 0-coverage samples
    if (verbose) loginfo(paste0("Filtering out regions with greater than ", filter2*100, "%% 0-coverage samples..."), logger = "LOGR")

    # GATK step 7 - filter out samples whose median coverage is outside the 2.5-97.5%ile range (SKIP)

    # GATK step 8 - Replace all zero coverages with target median
    if (verbose) loginfo("Setting all zero coverages to be the median over all samples...", logger = "LOGR")
    scaledrowmedian <- data[, rowMedians(as.matrix(.SD)), .SDcols = samples]
    data[, (samples) := lapply(.SD, function(values) {
        values[values == 0] <- scaledrowmedian[values == 0]
        values
    }), .SDcols = samples]

    # GATK step 9 - Winsorize
    if (verbose) loginfo(paste0("Clamping outliers outside the ", filter3*100, "%% and ", (1-filter3)*100, "%% percentiles to these values..."), logger = "LOGR")
    skip_winsorize = FALSE
    if (! skip_winsorize) {
        wins_range <- c(filter3, 1 - filter3)
        quantiles <- quantile(as.matrix(data), probs = wins_range)
        data[mask, (samples) := lapply(.SD, function(value) clip(value, quantiles[1], quantiles[2])), .SDcols = samples]
    }

    # GATK step 10 - Divide each column by column median
    if (verbose) loginfo("Dividing each sample by its median over all regions...", logger = "LOGR")
    data[mask, (samples) := lapply(.SD, function(values) values / median(values)), .SDcols = samples]

    # GATK step 11 - Take log2 of each coverage
    # Add any rows containing 0 or negative values to badrows mask, because can take log of these
    #badrows <- sort(union(badrows, which(as.matrix(data[, lapply(.SD, function(value) value <= 0), .SDcols = samples]), arr.ind = TRUE)[, 1]))
    if (verbose) loginfo("Converting to log2 domain...", logger = "LOGR")
    data[mask, (samples) := lapply(.SD, function(values) safelog2(values)), .SDcols = samples]

    # GATK step 12 - Calculate median of median of each sample, and subtract
    if (verbose) loginfo("Subtracting median of sample medians from each region...", logger = "LOGR")
    mom <- median(apply(data[mask], 2, median))
    data[mask, (samples) := lapply(.SD, function(values) values - mom), .SDcols = samples]

    # GATK step 13 - Calculate SVD
    if (verbose) loginfo("Computing SVD...", logger = "LOGR")
    .svd <- svd(data[mask])
    if (verbose) logdebug("Joliffe's criterion - %d", max(which(.svd$d > 0.7 * mean(.svd$d)))) # Joliffe's criterion
    if (verbose) loginfo("Completed", logger = "LOGR")
    list(mask = mask, u = .svd$u, data = data, rowmedian = rowmedian)
}


#' Do tangent normalisation on a tumour sample
#' @importFrom "logging" loginfo
#' @importFrom "data.table" setDT
#' @export
gatk_correct <- function(data, decomp, n_singular_values, verbose = TRUE) {
    setDT(data)
    samplenames <- colnames(data)
    data <- data[decomp$mask]

    if (n_singular_values > ncol(decomp$u)) {
        n_singular_values <- ncol(decomp$u)
    }

    if (n_singular_values < 1) {
        n_singular_values <- 1
    }

    if (verbose) loginfo("Converting to fractional read counts...", logger = "LOGR")
    data[, (samplenames) := lapply(.SD, function(values) values / sum(values)), .SDcols = samplenames]

    if (verbose) loginfo("Dividing regions by normal panel medians", logger = "LOGR")
    data[, (samplenames) := lapply(.SD, function(values) values / decomp$rowmedian[decomp$mask]), .SDcols = samplenames]

    if (verbose) loginfo("Converting to log2 domain...", logger = "LOGR")
    data[, (samplenames) := lapply(.SD, function(values) safelog2(values / median(values))), .SDcols = samplenames]

    if (verbose) loginfo("Standardizing values...", logger = "LOGR")
    t_standardized <- data#[, lapply(.SD, function(values) values - median(values)), .SDcols = samplenames]

    if (verbose) loginfo("Using SVD to remove noise...", logger = "LOGR")
    U <- decomp$u[, 1:n_singular_values]
    t_denoised <- t_standardized[, lapply(.SD, function(values) values - U %*% (t(U) %*% values)), .SDcols = samplenames]
    colnames(t_denoised) <- colnames(t_standardized)

    list(denoised = t_denoised, standardized = t_standardized)
}


#' GATK-style GC-correction of read count information
#' @importFrom "data.table" data.table
#' @export
gatk_gc_correct <- function(read_counts, gc_content, num_bins = 101, correlation_length = 0.02) {
    total <- sum(read_counts)
    which_bin <- round(gc_content * (num_bins - 1)) + 1
    correlation_decay_rate <- 1.0 / (correlation_length * num_bins)

    dt <- data.table(rc = read_counts, n = which_bin, gc = gc_content)
    dt <- dt[, .(med = median(rc), size = .N, gc = gc[1]), by = n]
    dt <- rbind(
        dt,
        data.table(n = setdiff(seq(1, num_bins), dt[, n]),
                   med = 1, size = 0, gc = setdiff(seq(1, num_bins), dt[, n]) / (num_bins - 1)),
        fill = TRUE
    )[order(n)]

    correction_factors <- sapply(dt[, n], function(bin) {
        weights <- dt[, size * exp(-abs(bin - n) * correlation_decay_rate)]
        sum(weights) / (weights %*% dt[, med])[1]
    })

    scaled <- read_counts * correction_factors[which_bin]
    scaled * total / sum(scaled)
}

