#' Finds the best fitting total copy number for the logR values in `logr_data_`,
#' given the `purity_` and `ploidy_` values.
#' @export
fit_total_copy_number <- function(logr_data_, purity_, ploidy_, max_cn_ = 10) {

    cn_states <- seq(0, max_cn_)
    # Calculate the expected logR value for each state
    ideal_logrs <- sapply(cn_states, function(n) {
        idealised_logr(n, 2, purity_, ploidy_)}
    )

    # Calculate the RMSE for each state
    logr_errors <- sapply(ideal_logrs, function(l) {
        sqrt(mean((logr_data_ - l)^2))
    })

    # Find the state which gives the smallest error
    min_index <- which.min(logr_errors)
    best_ntot <- test_ntots[min_index]

    return (
        list(total_copy_number=best_ntot,
             errors=logr_errors)
    )
}

#' Find the best fitting total copy number state for a vector of logR values,
#' for the given purity and ploidy values
#' @export
quick_fit_total_copy_number <- function(logr, purity, ploidy) {
    pointwise_total_cn <- calc_total_copynumber(logr, purity, ploidy)
    mean_total_cn <- mean(pointwise_total_cn)

    # define a test range of potential total copy number states as the mean ± 3 standard deviations
    test_range <- c(floor(mean_total_cn), ceiling(mean_total_cn)) + c(-3, 3) * sd(pointwise_total_cn)
    test_range[1] <- max(0, floor(test_range[1]))
    test_range[2] <- max(0, ceiling(test_range[2]))

    if (test_range[2] == 0) return (0)

    rmse <- function(v1, v2) sum(sqrt(mean((v1 - v2)^2)))

    test_states <- seq(test_range[1], test_range[2], 1)
    errors <- sapply(test_states, function(i) {
        rmse(idealised_logr(i, 2, purity, ploidy), logr)
    })

    return (test_states[which.min(errors)])
}


#' Finds the best fitting minor copy number for the VAF values in `vaf_data_`,
#' given the `purity_` and `total_cn` values.
#' @export
fit_minor_copy_number <- function(vaf_data_, purity_, total_cn) {
    # Find best fitting minor copy state
    max_nmin = as.integer(floor(total_cn/2))

    test_nmins <- seq(0, max_nmin)
    ideal_vafs <- sapply(test_nmins, function(n) {
        cnpipe::idealised_vaf(ta=(total_cn - n), tb=n, ha=1, hb=1, purity=purity_)
    })

    vaf_errors <- sapply(ideal_vafs, function(v) {
        sqrt(mean((mirror(vaf_data_) - v)^2))
    })

    min_index <- which.min(vaf_errors)
    best_nmin <- test_nmins[min_index]

    return (
        list(minor_copy_number=best_nmin,
             errors=vaf_errors)
    )
}
