#' Clip numeric vector so that its lowest element is >= low,
#' and its highest is <= high
clip <- function(v, low, high) {
    pmax(low, pmin(v, high))
}

safelog2 <- function(x) {
    clip(log2(x), -31, 31)
}

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
    loginfo(paste0("Data dimensions:", nrow(data), " regions; ", ncol(data), " samples.\n"), logger = "LOGR")
    loginfo("Making sure data is numeric...", logger = "LOGR")
    data[, (samples) := lapply(.SD, function(value) as.numeric(value)), .SDcols = samples]

    # Secret undocumented GATK step 1A - convert read counts to fractional read counts
    loginfo("Converting to fractional read counts...", logger = "LOGR")
    data[, (samples) := lapply(.SD, function(value) value / sum(value)), .SDcols = samples]

    # GATK step 2 - Calculate target medians
    loginfo("Calculating median across all samples...", logger = "LOGR")
    rowmedian = data[, rowMedians(as.matrix(.SD)), .SDcols = samples]

    # GATK step 3 - Apply 1st filter - remove intervals with median across all samples below 10%ile
    loginfo(paste0("Filtering out regions with a median across all samples below the ", filter1*100, "%% quantile..."), logger = "LOGR")
    mask <- (rowmedian >= quantile(rowmedian, filter1) & rowmedian != 0) # filters out bad rows - a vector of length == nrow(data),
    # TRUE for each row to retain, FALSE for each row to remove

    # GATK step 4: Divide each row by row median
    #rowmedian[rowmedian == 0] <- 1  # dummy count to avoid dividing by zero
    loginfo("Dividing each region by its median across all samples...", logger = "LOGR")
    data[mask] <- data[mask] / rowmedian[mask]

    # GATK step 5 - filter out samples with >5% 0-coverage targets (SKIP)

    # GATK step 6 - filter out targets with >5% 0-coverage samples
    loginfo(paste0("Filtering out regions with greater than ", filter2*100, "%% 0-coverage samples..."), logger = "LOGR")
    mask <- mask & (apply(data, 1, function(v) sum(v==0) / length(v)) <= filter2)

    # GATK step 7 - filter out samples whose median coverage is outside the 2.5-97.5%ile range (SKIP)

    # GATK step 8 - Replace all zero coverages with target median
    loginfo("Setting all zero coverages to be the median over all samples...", logger = "LOGR")
    scaledrowmedian <- data[, rowMedians(as.matrix(.SD)), .SDcols = samples]
    data[, (samples) := lapply(.SD, function(values) {
        values[values == 0] <- scaledrowmedian[values == 0]
        values
    }), .SDcols = samples]

    # GATK step 9 - Winsorize
    loginfo(paste0("Clamping outliers outside the ", filter3*100, "%% and ", (1-filter3)*100, "%% percentiles to these values..."), logger = "LOGR")
    skip_winsorize = FALSE
    if (! skip_winsorize) {
        wins_range <- c(filter3, 1 - filter3)
        quantiles <- quantile(as.matrix(data), probs = wins_range)
        data[mask, (samples) := lapply(.SD, function(value) clip(value, quantiles[1], quantiles[2])), .SDcols = samples]
    }

    # GATK step 10 - Divide each column by column median
    loginfo("Dividing each sample by its median over all regions...", logger = "LOGR")
    data[mask, (samples) := lapply(.SD, function(values) values / median(values)), .SDcols = samples]

    # GATK step 11 - Take log2 of each coverage
    # Add any rows containing 0 or negative values to badrows mask, because can take log of these
    #badrows <- sort(union(badrows, which(as.matrix(data[, lapply(.SD, function(value) value <= 0), .SDcols = samples]), arr.ind = TRUE)[, 1]))
    loginfo("Converting to log2 domain...", logger = "LOGR")
    data[mask, (samples) := lapply(.SD, function(values) safelog2(values)), .SDcols = samples]

    # GATK step 12 - Calculate median of median of each sample, and subtract
    loginfo("Subtracting median of sample medians from each region...", logger = "LOGR")
    mom <- median(apply(data[mask], 2, median))
    data[mask, (samples) := lapply(.SD, function(values) values - mom), .SDcols = samples]

    # GATK step 13 - Calculate SVD
    loginfo("Computing SVD...", logger = "LOGR")
    .svd <- svd(data[mask])
    #k <- max(which(.svd$d > 0.7 * mean(.svd$d))) # Joliffe's criterion
    k <- min(ncol(data), nsvs)
    u <- .svd$u[, 1:k]
    u_pinv <- MASS::ginv(u)
    loginfo("Completed", logger = "LOGR")
    list(mask = mask, u = u, data = data, rowmedian = rowmedian)
}


#' Do tangent normalisation on a tumour sample
#' @importFrom "logging" loginfo
#' @importFrom "data.table" setDT
#' @export
gatk_correct <- function(data, decomp, verbose = TRUE) {
    setDT(data)
    samplenames <- colnames(data)
    data <- data[decomp$mask]

    loginfo("Converting to fractional read counts...", logger = "LOGR")
    data[, (samplenames) := lapply(.SD, function(values) values / sum(values)), .SDcols = samplenames]

    loginfo("Dividing regions by normal panel medians", logger = "LOGR")
    data[, (samplenames) := lapply(.SD, function(values) values / decomp$rowmedian[decomp$mask]), .SDcols = samplenames]

    loginfo("Converting to log2 domain...", logger = "LOGR")
    data[, (samplenames) := lapply(.SD, function(values) safelog2(values / median(values))), .SDcols = samplenames]

    loginfo("Standardizing values...", logger = "LOGR")
    t_standardized <- data[, lapply(.SD, function(values) values - median(values)), .SDcols = samplenames]

    loginfo("Using SVD to remove noise...", logger = "LOGR")
    t_denoised <- t_standardized[, lapply(.SD, function(values) values - decomp$u %*% (t(decomp$u) %*% values)), .SDcols = samplenames]

    list(denoised = t_denoised, standardized = t_standardized)
}
