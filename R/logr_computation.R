#' Revised logr Calculation
#' 2 steps: 1 - get host median coverage data
#'          2 - compute logr based on median host coverage
#' @importFrom "matrixStats" rowMedians
#' @export
get_median_host <- function(m) {
    if (!is.matrix(m)) m <- as.matrix(m)
    rowMedians(sweep(m, 2, colSums(m), '/'))
}

#' Quick calculation of logR from host and tumour coverage data
#' Example:
#' # get reference median host information:
#' host_index = ... (index of host data columns)
#' median_host <- data[, get_median_host(.SD), .SDcols = host_index]
#'
#' # Compute logR for tumour columns
#' tumour_columns = ... (names of tumour columns [result written here])
#' tumour_index = ... (index of tumour coverage data columns)
#' stopifnot(length(tumour_columns) == length(tumour_index))
#' data[, (tumour_columns) := quick_logr(.SD, median_host),
#'      .SDcols = tumour_index]
#' @importFrom "matrixStats" colMedians
#' @export
quick_logr <- function(m, rowmed) {
    if (!is.matrix(m)) m <- as.matrix(m)
    stopifnot(length(rowmed) == nrow(m))
    m <- sweep(m, 2, colSums(m), '/')
    m <- sweep(m, 1, rowmed, '/')
    m[m==0] <- 1e-16
    m <- log2(m)
    cm <- colMedians(m)
    sweep(m, 2, cm, '-')
}
