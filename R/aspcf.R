#' If the penalty is set, then the segmentation is done to that penalty,
#' regardless of how many segments it ends up with
modded_aspcf <- function (ASCATobj, selectsamples = 1:length(ASCATobj$samples), 
          ascat.gg = NULL, penalty = 25, out.dir = ".", out.prefix = "") 
{
    madWins <- ASCAT:::madWins
    predictGermlineHomozygousStretches <- ASCAT:::predictGermlineHomozygousStretches
    fastAspcf <- ASCAT:::fastAspcf
    exactPcf <- ASCAT:::exactPcf
    fillNA <- ASCAT:::fillNA
    gg = NULL
    if (!is.null(ascat.gg)) {
        gg = ascat.gg$germlinegenotypes
    }
    else {
        gg = ASCATobj$Germline_BAF < 0.3 | ASCATobj$Germline_BAF > 
            0.7
    }
    ghs = predictGermlineHomozygousStretches(ASCATobj$chr, gg)
    segmentlengths = penalty
    segmentlengths = segmentlengths[segmentlengths >= penalty]
    Tumor_LogR_segmented = matrix(nrow = dim(ASCATobj$Tumor_LogR)[1], 
                                  ncol = dim(ASCATobj$Tumor_LogR)[2])
    rownames(Tumor_LogR_segmented) = rownames(ASCATobj$Tumor_LogR)
    colnames(Tumor_LogR_segmented) = colnames(ASCATobj$Tumor_LogR)
    Tumor_BAF_segmented = list()
    for (sample in selectsamples) {
        print.noquote(paste("Sample ", ASCATobj$samples[sample], 
                            " (", sample, "/", length(ASCATobj$samples), ")", 
                            sep = ""))
        logrfilename = file.path(out.dir, paste(out.prefix, 
                                                ASCATobj$samples[sample], ".LogR.PCFed.txt", sep = ""))
        baffilename = file.path(out.dir, paste(out.prefix, ASCATobj$samples[sample], 
                                               ".BAF.PCFed.txt", sep = ""))
        logRPCFed = numeric(0)
        bafPCFed = numeric(0)
        for (segmentlength in segmentlengths) {
            logRPCFed = numeric(0)
            bafPCFed = numeric(0)
            tbsam = ASCATobj$Tumor_BAF[, sample]
            names(tbsam) = rownames(ASCATobj$Tumor_BAF)
            homosam = gg[, sample]
            for (chrke in 1:length(ASCATobj$chr)) {
                lr = ASCATobj$Tumor_LogR[ASCATobj$chr[[chrke]], 
                                         sample]
                lrwins = vector(mode = "numeric", length = length(lr))
                lrwins[is.na(lr)] = NA
                lrwins[!is.na(lr)] = madWins(lr[!is.na(lr)], 
                                             2.5, 25)$ywin
                baf = tbsam[ASCATobj$chr[[chrke]]]
                homo = homosam[ASCATobj$chr[[chrke]]]
                Select_het <- !homo & !is.na(homo) & !is.na(baf) & 
                    !is.na(lr)
                bafsel = baf[Select_het]
                bafselwinsmirrored = madWins(ifelse(bafsel > 
                                                        0.5, bafsel, 1 - bafsel), 2.5, 25)$ywin
                bafselwins = ifelse(bafsel > 0.5, bafselwinsmirrored, 
                                    1 - bafselwinsmirrored)
                indices = which(Select_het)
                logRaveraged = NULL
                if (length(indices) != 0) {
                    averageIndices = c(1, (indices[1:(length(indices) - 
                                                          1)] + indices[2:length(indices)])/2, length(lr) + 
                                           0.01)
                    startindices = ceiling(averageIndices[1:(length(averageIndices) - 
                                                                 1)])
                    endindices = floor(averageIndices[2:length(averageIndices)] - 
                                           0.01)
                    if (length(indices) == 1) {
                        startindices = 1
                        endindices = length(lr)
                    }
                    nrIndices = endindices - startindices + 1
                    logRaveraged = vector(mode = "numeric", length = length(indices))
                    for (i in 1:length(indices)) {
                        if (is.na(endindices[i])) {
                            endindices[i] = startindices[i]
                        }
                        logRaveraged[i] = mean(lrwins[startindices[i]:endindices[i]], 
                                               na.rm = T)
                    }
                }
                if (length(logRaveraged) > 0) {
                    logRASPCF = NULL
                    bafASPCF = NULL
                    if (length(logRaveraged) < 6) {
                        logRASPCF = rep(mean(logRaveraged), length(logRaveraged))
                        bafASPCF = rep(mean(bafselwins), length(logRaveraged))
                    }
                    else {
                        PCFed = fastAspcf(logRaveraged, bafselwins, 
                                          6, segmentlength)
                        logRASPCF = PCFed$yhat1
                        bafASPCF = PCFed$yhat2
                    }
                    names(bafASPCF) = names(indices)
                    logRc = numeric(0)
                    for (probe in 1:length(logRASPCF)) {
                        if (probe == 1) {
                            logRc = rep(logRASPCF[probe], indices[probe])
                        }
                        if (probe == length(logRASPCF)) {
                            logRc = c(logRc, rep(logRASPCF[probe], 
                                                 length(lr) - indices[probe]))
                        }
                        else if (logRASPCF[probe] == logRASPCF[probe + 
                                                               1]) {
                            logRc = c(logRc, rep(logRASPCF[probe], 
                                                 indices[probe + 1] - indices[probe]))
                        }
                        else {
                            d = numeric(0)
                            totall = indices[probe + 1] - indices[probe]
                            for (bp in 0:(totall - 1)) {
                                dis = sum(abs(lr[(1:bp) + indices[probe]] - 
                                                  logRASPCF[probe]), na.rm = T)
                                if (bp != totall) {
                                    dis = sum(dis, sum(abs(lr[((bp + 1):totall) + 
                                                                  indices[probe]] - logRASPCF[probe + 
                                                                                                  1]), na.rm = T), na.rm = T)
                                }
                                d = c(d, dis)
                            }
                            breakpoint = which.min(d) - 1
                            logRc = c(logRc, rep(logRASPCF[probe], 
                                                 breakpoint), rep(logRASPCF[probe + 1], 
                                                                  totall - breakpoint))
                        }
                    }
                    logRd = numeric(0)
                    seg = rle(logRc)$lengths
                    startprobe = 1
                    endprobe = 0
                    for (i in 1:length(seg)) {
                        endprobe = endprobe + seg[i]
                        level = mean(lr[startprobe:endprobe], na.rm = T)
                        logRd = c(logRd, rep(level, seg[i]))
                        startprobe = startprobe + seg[i]
                    }
                    logRPCFed = c(logRPCFed, logRd)
                    bafPCFed = c(bafPCFed, bafASPCF)
                }
                else {
                    level = mean(lr, na.rm = T)
                    reps = length(lr)
                    logRPCFed = c(logRPCFed, rep(level, reps))
                }
                homsegs = ghs[[sample]][ghs[[sample]][, 1] == 
                                            chrke, ]
                startchr = min(ASCATobj$chr[[chrke]])
                endchr = max(ASCATobj$chr[[chrke]])
                if (length(homsegs) == 3) {
                    homsegs = t(as.matrix(homsegs))
                }
                if (!is.null(homsegs) && !is.na(homsegs) && 
                    dim(homsegs)[1] != 0) {
                    for (i in 1:dim(homsegs)[1]) {
                        startpos = max(homsegs[i, 2], startchr)
                        endpos = min(homsegs[i, 3], endchr)
                        startpos2 = max(homsegs[i, 2] - 100, startchr)
                        endpos2 = min(homsegs[i, 3] + 100, endchr)
                        startpos3 = max(homsegs[i, 2] - 5, startchr)
                        endpos3 = min(homsegs[i, 3] + 5, endchr)
                        towins = ASCATobj$Tumor_LogR[startpos2:endpos2, 
                                                     sample]
                        winsed = madWins(towins[!is.na(towins)], 
                                         2.5, 25)$ywin
                        pcfed = vector(mode = "numeric", length = length(towins))
                        pcfed[!is.na(towins)] = exactPcf(winsed, 
                                                         6, floor(segmentlength/4))
                        pcfed2 = pcfed[(startpos3 - startpos2 + 
                                            1):(endpos3 - startpos2 + 1)]
                        dif = abs(pcfed2 - logRPCFed[startpos3:endpos3])
                        if (!is.na(dif) && sum(dif > 0.3) > 5) {
                            logRPCFed[startpos3:endpos3] = ifelse(dif > 
                                                                      0.3, pcfed2, logRPCFed[startpos3:endpos3])
                        }
                    }
                }
            }
            logRPCFed = fillNA(logRPCFed, zeroIsNA = TRUE)
            seg = rle(logRPCFed)$lengths
            logRPCFed = numeric(0)
            startprobe = 1
            endprobe = 0
            prevlevel = 0
            for (i in 1:length(seg)) {
                endprobe = endprobe + seg[i]
                level = mean(ASCATobj$Tumor_LogR[startprobe:endprobe, 
                                                 sample], na.rm = T)
                if (is.nan(level)) {
                    level = prevlevel
                }
                else {
                    prevlevel = level
                }
                logRPCFed = c(logRPCFed, rep(level, seg[i]))
                startprobe = startprobe + seg[i]
            }
            names(logRPCFed) = rownames(ASCATobj$Tumor_LogR)
            if (length(unique(logRPCFed)) < 800) {
                break
            }
        }
        write.table(logRPCFed, logrfilename, sep = "\t", col.names = F)
        write.table(bafPCFed, baffilename, sep = "\t", col.names = F)
        bafPCFed = as.matrix(bafPCFed)
        Tumor_LogR_segmented[, sample] = logRPCFed
        Tumor_BAF_segmented[[sample]] = 1 - bafPCFed
    }
    ASCATobj = list(Tumor_LogR = ASCATobj$Tumor_LogR, Tumor_BAF = ASCATobj$Tumor_BAF, 
                    Tumor_LogR_segmented = Tumor_LogR_segmented, Tumor_BAF_segmented = Tumor_BAF_segmented, 
                    Germline_LogR = ASCATobj$Germline_LogR, Germline_BAF = ASCATobj$Germline_BAF, 
                    SNPpos = ASCATobj$SNPpos, ch = ASCATobj$ch, chr = ASCATobj$chr, 
                    chrs = ASCATobj$chrs, samples = colnames(ASCATobj$Tumor_LogR), 
                    gender = ASCATobj$gender, sexchromosomes = ASCATobj$sexchromosomes, 
                    failedarrays = ascat.gg$failedarrays)
    return(ASCATobj)
}
