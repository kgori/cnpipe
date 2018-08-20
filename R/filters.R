variant_depth_filter <- function(dt, column_names, min_variant_reads, min_vaf) {
    mask <- rep(FALSE, nrow(dt))
    for (col in column_names) {
        mask <- mask | dt[, eval(as.name(col)) >= min_variant_reads & eval(as.name(sub(".nr", ".vaf", col))) >= min_vaf]
    }
    mask[is.na(mask)] <- FALSE
    mask
}

