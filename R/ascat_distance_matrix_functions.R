#' Make a distance matrix directly from the segmentation data
#' (this is how runASCAT does it)
#' @param ascat_data ASCAT segmentation object
#' @export
make_distance_matrix <- function(ascat_data) {
    .make_segments <- function (r, b)
    {
        m = matrix(ncol = 2, nrow = length(b))
        m[, 1] = r
        m[, 2] = b
        m = as.matrix(na.omit(m))
        pcf_segments = matrix(ncol = 3, nrow = dim(m)[1])
        colnames(pcf_segments) = c("r", "b", "length")
        index = 0
        previousb = -1
        previousr = 1e+10
        for (i in 1:dim(m)[1]) {
            if (m[i, 2] != previousb || m[i, 1] != previousr) {
                index = index + 1
                count = 1
                pcf_segments[index, "r"] = m[i, 1]
                pcf_segments[index, "b"] = m[i, 2]
            }
            else {
                count = count + 1
            }
            pcf_segments[index, "length"] = count
            previousb = m[i, 2]
            previousr = m[i, 1]
        }
        pcf_segments = as.matrix(na.omit(pcf_segments))[, , drop = FALSE]
        return(pcf_segments)
    }

    .create_distance_matrix <- function (segments, gamma)
    {
        s = segments
        psi_pos = seq(1, 6, 0.05)
        rho_pos = seq(0.1, 1.05, 0.01)
        d = matrix(nrow = length(psi_pos), ncol = length(rho_pos))
        rownames(d) = psi_pos
        colnames(d) = rho_pos
        dmin = 1e+20
        for (i in 1:length(psi_pos)) {
            psi = psi_pos[i]
            for (j in 1:length(rho_pos)) {
                rho = rho_pos[j]
                nA = (rho - 1 - (s[, "b"] - 1) * 2^(s[, "r"]/gamma) *
                          ((1 - rho) * 2 + rho * psi))/rho
                nB = (rho - 1 + s[, "b"] * 2^(s[, "r"]/gamma) *
                          ((1 - rho) * 2 + rho * psi))/rho
                nMinor = NULL
                if (sum(nA, na.rm = T) < sum(nB, na.rm = T)) {
                    nMinor = nA
                }
                else {
                    nMinor = nB
                }
                d[i, j] = sum(abs(nMinor - pmax(round(nMinor), 0))^2 *
                                  s[, "length"] * ifelse(s[, "b"] == 0.5, 0.05,
                                                         1), na.rm = T)
            }
        }
        return(d)
    }

    dmlist <- lapply(1:dim(ascat_data$Tumor_LogR)[2], function(arraynr) {
        msg <- sprintf("%s - [ %d / %d ]\n", ascat_data$samples[arraynr], arraynr, length(ascat_data$samples))
        loginfo(msg)
        bafsegmented = ascat_data$Tumor_BAF_segmented[[arraynr]][, , drop = FALSE]
        names(bafsegmented) = rownames(ascat_data$Tumor_BAF_segmented[[arraynr]])
        lrrsegmented = ascat_data$Tumor_LogR_segmented[, arraynr]
        names(lrrsegmented) = rownames(ascat_data$Tumor_LogR_segmented)
        SNPpos <- ascat_data$SNPpos

        b = bafsegmented
        r = lrrsegmented[names(bafsegmented)]
        SNPposhet = SNPpos[names(bafsegmented), ]
        autoprobes = !(SNPposhet[, 1] %in% ascat_data$sexchromosomes)
        b2 = b[autoprobes]
        r2 = r[autoprobes]

        s <- .make_segments(r2, b2)

        dm <- .create_distance_matrix(s, 1.0) # gamma = 1.0 for NGS
        return (dm)
    })
    names(dmlist) <- ascat_data$samples
    dmlist
} # make_distance_matrix


#' Find peaks in the ASCAT sunrise plot data.
#' @param dm = Distance matrix from ascat_result$distance_matrix
#' @param nbr_size = Neighbourhood size for peakfinding (e.g. 5x5 matrix centred on peak)
# @importFrom "raster" raster extent focal
#' @export
# peakfind <- function(dm, nbr_size=5) {
#     ## Convert distance matrix to a raster object
#     # The DM is transposed and reversed along 1 axis to make sure the
#     # ploidy and purity axes line up with the data in the raster object
#     r <- raster(t(log(dm))[rev(1:ncol(dm)), 1:nrow(dm)])
#     r@extent <- extent(c(range(as.numeric(rownames(dm))),
#                           range(as.numeric(colnames(dm)))))
#
#     ## Find the minimum value within the 25-cell neighborhood of each cell
#     f <- function(X) min(X, na.rm=TRUE)
#     ww <- matrix(1, nrow = nbr_size, ncol = nbr_size) ## Weight matrix for cells in moving window
#     localmin <- focal(r, fun = f, w = ww, pad = TRUE, padValue = NA)
#
#     ## Does each cell have the maximum value in its neighborhood?
#     r2 <- r == localmin
#
#     ## Get x-y coordinates of those cells that are local maxima
#     maxXY <- xyFromCell(r2, Which(r2==1, cells=TRUE))
#     colnames(maxXY) <- c("psi", "rho")
#     maxXY <- maxXY[maxXY[, 2] <= 1, , drop = FALSE] # filter out any purity estimates > 1
#     maxXY <- cbind(maxXY, apply(maxXY, 1, function(r) interpolate_peak(dm, r[1], r[2])))
#     colnames(maxXY)[3] <- "fit"
#     maxXY <- maxXY[order(maxXY[, 3]), , drop = FALSE] # sort by goodness of fit
#
#     list(maxXY, r, r2)
# }

find_box_index <- function(dm, psi, rho) {
    x0 <- findInterval(psi, as.numeric(rownames(dm)))
    x1 <- x0 + 1
    y0 <- findInterval(rho, as.numeric(colnames(dm)))
    y1 <- y0 + 1

    # m = Proportion "psi" is above x0, in the interval [psi_x0, psi_x1)
    # n = Proportion "rho" is above y0, in the interval [rho_y0, rho_y1)
    psi_x0 <- as.numeric(rownames(dm)[x0])
    psi_x1 <- as.numeric(rownames(dm)[x1])
    rho_y0 <- as.numeric(colnames(dm)[y0])
    rho_y1 <- as.numeric(colnames(dm)[y1])
    m <- (psi - psi_x0) / (psi_x1 - psi_x0)
    n <- (rho - rho_y0) / (rho_y1 - rho_y0)
    list(x0=x0, x1=x1, y0=y0, y1=y1, m=m, n=n)
}

get_box_values <- function(dm, box) {
    a <- dm[box$x0, box$y0]
    b <- dm[box$x1, box$y0]
    c <- dm[box$x0, box$y1]
    d <- dm[box$x1, box$y1]
    list(a=a, b=b, c=c, d=d)
}

box_interpolate <- function(boxix, boxvals) {
    m <- boxix$m
    n <- boxix$n
    u <- boxvals$a * (1-m) + m * boxvals$b
    v <- boxvals$c * (1-m) + m * boxvals$d
    u * (1-n) + n*v
}

#' Interpolates the value of the ASCAT distance matrix at (psi, rho)
#' @param dm Matrix Distance matrix produced by ASCAT, or calculated from
#' an ASCAT segmentation object using the function \code{make_distance_matrix}
#' @export
interpolate_peak <- function(dm, psi, rho) {
    ix <- find_box_index(dm, psi, rho)
    vals <- get_box_values(dm, ix)
    box_interpolate(ix, vals)
}

#' @export
plot_dm <- function(dm, peaks = NULL) {
    ASCAT::ascat.plotSunrise(dm, 0, 0)
    if (!is.null(peaks)) {
        ploidy_min <- as.numeric(rownames(dm)[1])
        ploidy_max <- as.numeric(rownames(dm)[nrow(dm)])
        purity_min <- as.numeric(colnames(dm)[1])
        purity_max <- as.numeric(colnames(dm)[ncol(dm)])
        cex.crosses <- (function(v) (2/v)/max(1/v))(peaks[,3])
        points((peaks[,1] - ploidy_min)/(ploidy_max - 1),
               (peaks[,2] - purity_min)/(1/purity_max),
               col = "green", pch = 4,
               cex = cex.crosses)
    }
}

#' Extract a single sample from an ascat data object
#' @export
ascat_single_sample <- function(ascat_data, samplename) {
    stopifnot(samplename %in% ascat_data$samples)
    index <- which(ascat_data$samples == samplename)
    new <- list()
    new$Germline_BAF         <- ascat_data$Germline_BAF[, index, drop = FALSE]
    new$Germline_LogR        <- ascat_data$Germline_LogR[, index, drop = FALSE]
    new$Tumor_BAF            <- ascat_data$Tumor_BAF[, index, drop = FALSE]
    new$Tumor_LogR           <- ascat_data$Tumor_LogR[, index, drop = FALSE]
    new$Tumor_BAF_segmented  <- list(ascat_data$Tumor_BAF_segmented[[index]])
    new$Tumor_LogR_segmented <- ascat_data$Tumor_LogR_segmented[, index, drop = FALSE]
    new$gender               <- ascat_data$gender[index]
    new$samples              <- ascat_data$samples[index]
    new$SNPpos               <- ascat_data$SNPpos
    new$ch                   <- ascat_data$ch
    new$chr                  <- ascat_data$chr
    new$chrs                 <- ascat_data$chrs
    new$sexchromosomes       <- ascat_data$sexchromosomes
    new$failedarrays         <- ascat_data$failedarrays
    new
}
