#' Load the two tumour files that ASCAT worked on
#' @importFrom "GenomicRanges" GRanges
#' @importFrom "data.table" fread
#' @export
karyoplot_load_baflogr <- function(filename.baf, filename.logr) {
    baf <- fread(filename.baf, sep='\t', drop = 1)
    logr <- fread(filename.logr, sep='\t', drop = 1)
    colnames(baf)[3] <- "BAF"
    colnames(logr)[3] <- "LOGR"
    baf[, c("start", "end") := pos]
    baf[, "LOGR" := logr$LOGR]
    baf[, pos := NULL]
    GRanges(baf)
}

#' Load the segmentation that ASCAT produced
#' @importFrom "GenomicRanges" GRanges
#' @importFrom "data.table" fread
#' @export
karyoplot_load_segments <- function(filename) {
    cn.segments <- fread(filename, sep='\t')
    colnames(cn.segments)[3:4] <- c("start", "end")
    GRanges(cn.segments)
}

#' Build a GRanges object describing the dog genome
#' @importFrom "regioneR" toGRanges
#' @export
karyoplot_dog_genome <- function() {
    toGRanges(
        data.frame(
            chr = c(
                "1", "2", "3", "4", "5",
                "6", "7", "8", "9", "10",
                "11", "12", "13", "14", "15",
                "16", "17", "18", "19", "20",
                "21", "22", "23", "24", "25",
                "26", "27", "28", "29", "30",
                "31", "32", "33", "34", "35",
                "36", "37", "38", "X"
            ),
            start = rep(1, 39),
            end = c(
                122678785, 85426708, 91889043, 88276631, 88915250,
                77573801, 80974532, 74330416, 61074082, 69331447,
                74389097, 72498081, 63241923, 60966679, 64190966,
                59632846, 64289059, 55844845, 53741614, 58134056,
                50858623, 61439934, 52294480, 47698779, 51628933,
                38964690, 45876710, 41182112, 41845238, 40214260,
                39895921, 38810281, 31377067, 42124431, 26524999,
                30810995, 30902991, 23914537, 123869142
            )
        )
    )
}

#' Plot the karyplot diagram
#' @importFrom "karyoploteR" getDefaultPlotParams plotKaryotype kpAddCytobandsAsLine
#'                           kpAddChromosomeNames kpAddBaseNumbers kpAddMainTitle
#'                           kpAxis kpPoints kpSegments kpAbline
#' @export
karyoplot <- function(baflogr.data, segments.data, genome, chromosomes, title) {
    par(lend=1)
    chrom <- chromosomes
    pt <- 1
    pp <- getDefaultPlotParams(plot.type = pt)
    # pp$data1inmargin <- 5
    pp$ideogramlateralmargin <- 0.05 # default: 0.003
    kp <- plotKaryotype(genome = genome, chromosomes = chrom,
                        plot.type=pt, ideogram.plotter = NULL,
                        labels.plotter = NULL, plot.params = pp)
    kpAddCytobandsAsLine(kp)
    kpAddChromosomeNames(kp, srt=45)
    kpAddBaseNumbers(kp)
    kpAddMainTitle(kp, main = title, cex=1.2)

    # LOGR
    logr.r0 <- 0.0
    logr.r1 <- 0.25
    logr.ymin <- -3
    logr.ymax <- 5
    kpAxis(kp, r0=logr.r0, r1=logr.r1, tick.pos = c(logr.ymin, 0, logr.ymax), ymin=logr.ymin, ymax=logr.ymax)
    kpPoints(kp, data=baflogr.data, y=baflogr.data$LOGR, cex = 0.2, col="black", r0=logr.r0, r1=logr.r1, ymin=logr.ymin, ymax=logr.ymax)

    # BAF
    baf.r0 <- 0.75
    baf.r1 <- 1.0
    kpAxis(kp, r0=baf.r0, r1=baf.r1)
    kpPoints(kp, data=baflogr.data, y=baflogr.data$BAF, cex = 0.2, col="black", r0=baf.r0, r1=baf.r1)

    # Plot copy number segments in their own section
    cn.r0 = 0.3
    cn.r1 = 0.7
    cn.ymin <- 0
    cn.ymax <- 6
    kpAxis(kp, r0=cn.r0, r1=cn.r1, ymin=cn.ymin, ymax=cn.ymax, tick.pos = seq(cn.ymin, cn.ymax), side=2, gr)
    kpSegments(kp, data=segments.data, y0=(segments.data$nMajor+segments.data$nMinor)+0.05, y1=(segments.data$nMajor+segments.data$nMinor)+0.05, ymin=cn.ymin, ymax=cn.ymax, lwd=2, col = "red", r0 = cn.r0, r1 = cn.r1)
    kpSegments(kp, data=segments.data, y0=(segments.data$nMinor)-0.05, y1=(segments.data$nMinor)-0.05, ymin=cn.ymin, ymax=cn.ymax, lwd=2, col = "green", r0 = cn.r0, r1 = cn.r1)
    kpSegments(kp, data=segments.data, y0=(segments.data$nMajor), y1=(segments.data$nMajor), ymin=cn.ymin, ymax=cn.ymax, lwd=2, col = "blue", r0 = cn.r0, r1 = cn.r1)

    for (c in chrom) {
        breakpoints <- segments.data[segments.data@seqnames==c]@ranges@start
        breakpoints <- breakpoints[2:length(breakpoints)]
        kpAbline(kp, chr = c, v=breakpoints, col="gray", lty = 1, lwd = .5)
    }
    kp
}
