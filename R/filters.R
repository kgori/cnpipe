#' @export
variant_depth_filter <- function(dt, variant_column_names, min_variant_reads, min_vaf) {
    mask <- rep(FALSE, nrow(dt))
    if (length(min_variant_reads) == length(variant_column_names)) {
        thresholds <- data.frame(name = variant_column_names,
                                 threshold = min_variant_reads,
                                 stringsAsFactors = FALSE)
    } else {
        thresholds <- data.frame(name = variant_column_names,
                                 threshold = rep(min_variant_reads[1], length(variant_column_names)),
                                 stringsAsFactors = FALSE)
    }
    for (i in 1:nrow(thresholds)) {
        nv <- thresholds[i, ][["name"]]
        vaf <- sub(".nv", ".vaf", nv)
        min_nv <- thresholds[i, ][["threshold"]]
        mask <- mask | dt[, eval(as.name(nv)) >= min_nv | eval(as.name(vaf)) >= min_vaf]
    }
    mask[is.na(mask)] <- FALSE
    mask
}
