#' Implementation of PCF algorithm from
#' Copynumber: Efficient algorithms for
#' single and multi-track copy number
#' segmentation
#' (not very efficient version)
#' example:
#' data = c(rnorm(495, 0, 0.1), rnorm(10, 0.5, 0.1), rnorm(495, 0, 0.1), rnorm(100, 1, 0.1), rnorm(1000, 1, 0.1))
#' pcf(data)
pcf <- function(y, gamma = 40) {
    p = length(y)
    a0 = numeric(0)
    e0 = 0
    aprev = a0
    eprev = e0
    gamma = gamma * var(y)
    t <- rep(0, p) # backtrack
    logdebug(paste0("a0=",a0," e0=",e0," aprev=",aprev," eprev=", eprev," gamma=",gamma))
    for (k in seq(1, p)) {
        anow <- c(aprev, 0) + y[k]
        dnow <- anow^2 / (1:k - k - 1)
        score <- (dnow + eprev + gamma)
        min_index <- which.min(score)
        enow <- c(eprev, score[min_index])
        t[k] <- min_index
        loginfo(paste0("anow=",anow," enow=",enow," dnow=", dnow," score=",score))
        aprev <- anow
        eprev <- enow
    }

    segstarts <- unique(t)
    segends <- c(segstarts[2:length(segstarts)]-1, p)
    data.frame(startpos=segstarts, endpos=segends)
}

# ys is a matrix of datapoints (rows) x samples (cols)
# y <- cbind(c(rnorm(1000, 0, 0.1), rnorm(100, 1, 0.1), rnorm(1000, 0, 0.1)), c(rnorm(1000, 0, 0.1), rnorm(100, 1, 0.1), rnorm(400, 0, 0.1), rnorm(100, 1, 0.1), rnorm(500, 0, 0.1)))
multipcf <- function(ys, gamma = 40) {
    p = nrow(ys)
    a0 = matrix(0, nrow = nrow(ys), ncol = ncol(ys))
    e0 = 0
    aprev = a0
    aprev[1,] <- ys[1,]
    anow = a0
    eprev = e0
    t <- rep(0, p) # backtrack
    for (k in seq(1, p)) {
        anow[1:(k+1)] <- aprev[1:(k+1), ] + ys[k,]
    }
}


#' Applies a prefilter to a PCF segmentation
#' @param dta The data that was segmented
#' @param seg The result of segmenting dta
#' @return A filtered segmentation
prefilter <- function(dta, seg, pthresh=0.3, med_diff=0.5) {
    stopifnot(pthresh >= 0)
    stopifnot(pthresh <= 1)
    stopifnot(med_diff %in% c(0.5, 1.0))
    stopifnot("total_cn_tetraploid" %in% colnames(dta))
    stopifnot("START" %in% colnames(dta))

    # Annotate the table (dta) with the segment ID (bk)
    dta[, bk := rep(seq(1, seg$nIntervals), times = seg$Lengde)]
    dta[, segment_median := median(updated_cn), by = bk]

    # Collect some stats for each segment, and its successor
    test_table <- dta[, .(START=min(START), mn=mean(updated_cn), sqsem=var(updated_cn) / .N, imdn=round(median(updated_cn)), mdn=0.5*round(2*median(updated_cn)), N=.N), by = bk]
    test_table[, c("bk_next", "mn_next", "sqsem_next", "imdn_next", "mdn_next", "N_next") := shift(.(bk, mn, sqsem, imdn, mdn, N), type="lead")]

    # Probability that the difference of the means is <0.5 by chance
    test_table[, pval := pnorm(0.5, abs(mn_next - mn), sqrt(sqsem + sqsem_next))]

    # Successive segment medians are different
    test_table[, medians_differ := mdn != mdn_next]
    test_table[, integer_medians_differ := imdn != imdn_next]

    approved_breaks <- if (med_diff == 0.5) {
        test_table[pval < pthresh & medians_differ == TRUE, bk_next]
    } else {
        test_table[pval < pthresh & integer_medians_differ == TRUE, bk_next]
    }

    is_breakpoint_start = dta[bk %in% approved_breaks, .(START=min(START)), by = bk]
    dta[, tmp := 0]
    dta[is_breakpoint_start, tmp := 1, on = c("START")]
    dta[, filtered_seg_id := cumsum(tmp) + 1]
    dta[, tmpI := .I]

    result <- dta[, .(Lengde=.N, sta=min(tmpI), mean=mean(total_cn_tetraploid)), by = filtered_seg_id]
    result <- list(Lengde = result[, Lengde],
                   sta = result[, sta],
                   mean = result[, mean],
                   nIntervals = dta[, max(filtered_seg_id)])

    dta[, tmp := NULL]
    dta[, tmpI := NULL]
    dta[, unfiltered_seg_id := bk]
    dta[, bk := NULL]
    dta[, segment_median := NULL]

    return (result)
}
