# Pack data into an ASCAT object

#' Helper function to get the sample name from a data.table
#' of tumour data
extract_sample_name <- function(dt) {
    dt[1, "samplename"][[1]]
}

#' Helper function to convert a tall data.table of genotyped breakpoints
#' into segment intervals. ('Tall' means the data for each sample
#' is stacked vertically (melted), as opposed to horizontally, in
#' separate columns (dcasted).
#' @param breaks - data.table with columns 'chr', 'pos', 'sample' and 'PASS',
#' where (chr,pos) are the coordinates of a breakpoint, sample is the name
#' of the sample the breakpoint appears in, and 'PASS' is a boolean indicating
#' whether the breakpoint exists in the sample
#' @param snpdata - data.table with columns 'chr', 'pos', where (chr,pos) are
#' SNP coordinates
#' @return - data.table of (chr,startpos,endpos) intervals, where startpos
#' and endpos are tied to available SNP positions from snpdata
convert_genotyped_breaks_to_intervals <- function(breaks, snpdata, allow_pass_in_any_sample = FALSE) {
    # Annotate each row of snpdata with the position of its preceding
    # SNP - used later to define interval end points
    snpdata[, prev := shift(pos), by = chr]

    # Collect all passing breakpoints
    setnames(breaks, old = "sample", new = "samplename", skip_absent = TRUE)
    if (allow_pass_in_any_sample) {
        passing <- breaks[breaks[, .(AllFail = all(PASS==FALSE)), by = .(chr, pos)][!(AllFail)], on = .(chr, pos)][, AllFail := NULL]
        passing <- passing[, .(chr, pos, samplename)]
    } else {
        passing <- breaks[(PASS), .(chr, pos, samplename)]
    }

    # Append first and last site of each chr as de facto breakpoints.
    # The dummy column, k, makes the upcoming join work
    first <- snpdata[, .SD[1], by = chr][, .(chr, pos, k=1)]
    last <- snpdata[, .SD[.N], by = chr][, .(chr, pos, k=1)]
    u <- passing[, .(samplename=unique(samplename), k=1)]
    setkey(first, k)
    setkey(last, k)
    setkey(u, k)
    first <- first[u, allow.cartesian=TRUE][,k:=NULL]
    last <- last[u, allow.cartesian=TRUE][,k:=NULL]
    passing <- rbindlist(list(passing, first, last))
    setorder(passing, samplename, chr, pos)

    # Convert breakpoints to contiguous intervals (ivls). a and b are temporaries
    # to make the sequence of joins clearer.
    # 1: Annotate with prev, the coordinate of the preceding SNP in snpdata (relative to pos)
    a = passing[snpdata, on = .(chr, pos), nomatch=0L]

    # 2: For each row, begin a segment at the current row's 'pos', and use the 'prev' of the *next*
    # row to define the end point of the segment
    b = a[, .(chr, startpos = pos, endpos = shift(prev, type = "lead")), by = samplename]

    # 3: Remove any segments that begin at the last SNP of the chromosome, and adjust preceding row to include final SNP
    b[is.na(shift(endpos, type="lead")) | is.na(endpos),
      endpos := ifelse(!is.na(endpos), shift(startpos, type = "lead"), NA)] # adjustment
    ivls = b[!is.na(endpos)] # removal

    return (ivls)
}



#' Output: data.frame [Nsites * Nsamples]
pack_tumour_logr <- function(sample_data_list) {
    df <- data.frame(sample_data_list[[1]]$T_logR)
    colnames(df) <- extract_sample_name(sample_data_list[[1]])
    if (length(sample_data_list) > 1) {
        for (i in 2:length(sample_data_list)) {
            sample_data <- sample_data_list[[i]]
            df <- cbind(df, sample_data$T_logR)
            colnames(df)[i] <- extract_sample_name(sample_data)
        }
    }
    df
}

#' Output: data.frame [Nsites * Nsamples]
pack_tumour_baf <- function(sample_data_list) {
    df <- data.frame(sample_data_list[[1]][, "T_corrected_vaf"])
    colnames(df) <- extract_sample_name(sample_data_list[[1]])
    if (length(sample_data_list) > 1) {
        for (i in 2:length(sample_data_list)) {
            sample_data <- sample_data_list[[i]]
            df <- cbind(df, sample_data[, "T_corrected_vaf"])
            colnames(df)[i] <- extract_sample_name(sample_data)
        }
    }
    df
}

#' Output: matrix [Nsites * Nsamples]
pack_tumour_logr_segmented <- function(sample_data_list, segmentation) {
    ret <- list()
    setkey(segmentation, samplename, chr, startpos, endpos)
    for (sample_data in sample_data_list) {
        samplename <- extract_sample_name(sample_data)
        sample_data[, c("chr", "startpos", "endpos") := .(CHROM, POS, POS)]
        setkeyv(sample_data, key(segmentation))
        logr_means <- foverlaps(sample_data, segmentation)[, mean_logR := mean(T_logR), by = .(CHROM, startpos, endpos)][, mean_logR]
        ret[[samplename]] <- logr_means
    }
    # Don't coerce colnames
    return (as.data.frame(ret, optional = TRUE))
}

#' Output: list [length=Nsamples] of matrix [NHetSites * 1]
#' rownames of matrix must correspond to the indices of SNPpos
#' that point to heterozygous sites
pack_tumour_baf_segmented <- function(sample_data_list, segmentation) {
    ret <- list()
    setkey(segmentation, samplename, chr, startpos, endpos)
    for (sample_data in sample_data_list) {
        # Here use the ASCAT heuristic for solving double-banding, with an adjustment:
        # the variance is calculated for the un-mirrored VAF data
        sample_data[, c("chr", "startpos", "endpos") := .(CHROM, POS, POS)]
        samplename <- extract_sample_name(sample_data)
        setkeyv(sample_data, key(segmentation))
        ol <- foverlaps(sample_data, segmentation)

        sd2 <- ol[, .(sd2 = getMad(T_corrected_vaf2)), by = .(chr, startpos, endpos)][, sqrt(mean(sd2^2))]
        mutab <- ol[, .(mu = mean(abs(mirror(T_corrected_vaf2) - 0.5))), by = .(chr, startpos, endpos)]
        mutab[sqrt(sd2^2+mu^2) < 2 * sd2, mu := 0]

        values <- ol[mutab, mean_vaf := mirror(mu + 0.5), on = .(chr, startpos, endpos)][, mean_vaf]
        ret[[length(ret)+1]] <- matrix(values, ncol=1, dimnames = list(seq_along(values)))[-1, 1, drop=FALSE]
    }
    return (ret)
}

#' Output: data.frame [Nsites * Nsamples]
pack_germline_logr <- function(sample_data_list) {
    ret <- as.data.frame(lapply(sample_data_list, function(sample_data) {
        rep(0, sample_data[,.N])
    }))
    names(ret) <- paste0(sapply(sample_data_list, extract_sample_name), ".host")
    ret
}

#' Output: data.frame [Nsites * Nsamples]
pack_germline_baf <- function(sample_data_list) {
    ret <- as.data.frame(lapply(sample_data_list, function(sample_data) {
        vals <- rep(0.5, sample_data[,.N])
        vals[1] <- 0
        vals
    }))
    names(ret) <- paste0(sapply(sample_data_list, extract_sample_name), ".host")
    ret
}

#' Output: data.frame [NSites * 2]
#' colnames are "CHROM" and "POS"
pack_snppos <- function(data) {
    df <- data.frame(
        CHROM = data$CHROM,
        POS = data$POS)
    df
}

#' Output: 2 * list [length=Nchrom], where list[[K]]
#' is the sequence 1:N, where N is the number of datapoints
#' coming from chromosome K
#' obj$ch and obj$chr are identical
#' obj$chrs is a character vector of each K
pack_chrs <- function(data) {
    ch <- list()
    data[, I := 1:.N, by = CHROM]
    chrs <- data[, unique(CHROM)]
    for (chrom in chrs) {
        ch[[length(ch) + 1]] <- data[CHROM==chrom, I]
    }
    return (list(ch=ch, chr=ch, chrs=chrs))
}

#' Output character vector [length=Nsamples]
pack_samples <- function(sample_data_list) {
    sapply(sample_data_list, extract_sample_name)
}

#' Output: character vector [length=Nsamples]
#' Sex of each sample ('XX' or 'XY' if sexchromosomes are 'X', 'Y')
pack_gender <- function(sample_data_list) {
    sapply(sample_data_list, function(nothing) {
        "XX"
    })
}

#' Output character vector [length=2]
pack_sexchromosomes <- function() {
    c('X', 'Y')
}


#' An ASCAT segmentation object is a list of length 14 containing:
#' 1  $Tumor_LogR
#' 2  $Tumor_BAF
#' 3  $Tumor_LogR_segmented
#' 4  $Tumor_BAF_segmented
#' 5  $Germline_LogR
#' 6  $Germline_BAF
#' 7  $SNPpos
#' 8  $ch
#' 9  $chr
#' 10 $chrs
#' 11 $samples
#' 12 $gender
#' 13 $sexchromosomes
#' 14 $failedarrays
#' @export
pack_ascat <- function(data, segmentation) {
    obj <- list()
    obj$Tumor_LogR <- pack_tumour_logr(data)
    obj$Tumor_BAF <- pack_tumour_baf(data)
    obj$Tumor_LogR_segmented <- pack_tumour_logr_segmented(data, segmentation)
    obj$Tumor_BAF_segmented <- pack_tumour_baf_segmented(data, segmentation)
    obj$Germline_LogR <- pack_germline_logr(data)
    obj$Germline_BAF <- pack_germline_baf(data)
    obj$SNPpos <- pack_snppos(data[[1]])
    ch <- pack_chrs(data[[1]])
    obj$ch <- ch$ch
    obj$chr <- ch$chr
    obj$chrs <- ch$chrs
    obj$samples <- pack_samples(data)
    obj$gender <- pack_gender(data)
    obj$sexchromosomes <- pack_sexchromosomes()
    obj$failedarrays <- NULL
    obj
}

pack_tumour_baf_segmented_plots <- function(sample_data_list, segmentation) {
    setkey(segmentation, samplename, chr, startpos, endpos)
    for (sample_data in sample_data_list) {
        pdf("~/Documents/projects/copynumber/pipeline_work_in_progress/CTVT/Oct_2019/16_10_19/double_banding_99siglevel_withlogR.pdf",
            width = 40, height = 12)
        par(mfrow = c(3, 1), mar = c(4,4,4,2))
        sample_data[, c("chr", "startpos", "endpos") := .(CHROM, POS, POS)]
        setkeyv(sample_data, key(segmentation))
        ol <- foverlaps(sample_data, segmentation)

        # Here use the ASCAT heuristic for solving double-banding
        sd2 <- ol[, .(sd2 = getMad(mirror(T_corrected_vaf2))), by = .(chr, startpos, endpos)][, sqrt(mean(sd2^2))]
        mutab <- ol[, .(mu = mean(abs(mirror(T_corrected_vaf2) - 0.5))), by = .(chr, startpos, endpos)]
        mutab[sqrt(sd2^2+mu^2) < 2 * sd2, mu := 0]

        # debug - plot
        ol[mutab, mu := mirror(mu + 0.5), on = .(chr, startpos, endpos)]
        segs = unique(ol[, .(startI = min(I), endI = max(I), mu=mean(T_logR)), by = .(chr, startpos, endpos)])
        plot(ol[chr=="1", .(I,T_logR)], pch = 20, col = scales::alpha("skyblue", 0.5), cex = .5, main="logR")
        segments(segs[chr=="1", startI], segs[chr=="1", mu], segs[chr=="1", endI], col = "blue", lwd = 4, lend = "butt")
        #segments(segs[chr=="1", startI], segs[chr=="1", 1-mu], segs[chr=="1", endI], col = "blue", lwd = 4)

        # Here use the ASCAT heuristic for solving double-banding
        sample_data[, c("chr", "startpos", "endpos") := .(CHROM, POS, POS)]
        setkeyv(sample_data, key(segmentation))
        ol <- foverlaps(sample_data, segmentation)

        sd2 <- ol[, .(sd2 = getMad(T_corrected_vaf2)), by = .(chr, startpos, endpos)][, sqrt(mean(sd2^2))]
        mutab <- ol[, .(mu = mean(abs(mirror(T_corrected_vaf2) - 0.5))), by = .(chr, startpos, endpos)]
        mutab[sqrt(sd2^2+mu^2) < 2 * sd2, mu := 0]

        # debug - plot
        ol[mutab, mu := mirror(mu + 0.5), on = .(chr, startpos, endpos)]
        segs = unique(ol[, .(startI = min(I), endI = max(I), mu), by = .(chr, startpos, endpos)])
        plot(ol[chr=="1", .(I, T_corrected_vaf2)], pch = 20, col = scales::alpha("red", 0.5), cex = .5, main="after - quick fix")
        segments(segs[chr=="1", startI], segs[chr=="1", mu], segs[chr=="1", endI], col = "green", lwd = 4, lend = "butt")
        segments(segs[chr=="1", startI], segs[chr=="1", 1-mu], segs[chr=="1", endI], col = "green", lwd = 4, lend = "butt")

        # my method
        sample_data[, c("chr", "startpos", "endpos") := .(CHROM, POS, POS)]
        setkeyv(sample_data, key(segmentation))
        ol <- foverlaps(sample_data, segmentation)
        rm(segs)

        tester <- ol[chr=="1"]
        tester[, mu:= dotest(T_corrected_vaf2), by = .(chr, startpos, endpos)]
        plot(tester[chr=="1", .(I, T_corrected_vaf2)], pch = 20, col = scales::alpha("red", 0.5), cex = .5, main="after - statistical fix")
        segs = unique(tester[, .(startI = min(I), endI = max(I), mu), by = .(chr, startpos, endpos)])
        segments(segs[chr=="1", startI], segs[chr=="1", mu], segs[chr=="1", endI], col = "green", lwd = 4, lend = "butt")
        segments(segs[chr=="1", startI], segs[chr=="1", 1-mu], segs[chr=="1", endI], col = "green", lwd = 3, lend = "butt")
        dev.off()
    }
}
