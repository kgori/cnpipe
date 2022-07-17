#' Estimate the host and tumour components from a host-contaminated tumour sample
#' @param tumour_ref number of REF reads observed in host-contaminated tumour, at current site
#' @param tumour_alt number of ALT reads observed in host-contaminated tumour, at current site
#' @param host_ref number of REF reads observed in matched host, at current site
#' @param host_alt number of ALT reads observed in matched host, at current site
#' @param logr log of ratio between tumour coverage and matched host coverage, at current site
#' @param purity Estimate of tumour purity (aka aberrant cell count)
#' @param plot Should a plot be drawn?
#' @importFrom "mixtools" ellipse
#' @return A list with the following elements:
#' * par = optimised parameters (\eqn{p(H | A)}, \eqn{p(H | B)})
#' * cov = estimated covariance of optimised parameters
#'       (Fisher information - inverse of Hessian)
#' * counts = Matrix of inferred read counts, divided by Host/Tumour, A/B-allele
#' * pure_tumour_vaf = Estimate of tumour VAF in absence of host
#' * pure_host_vaf = Estimate of host VAF in absence of tumour
#' * sample_vaf_with_het_host = An estimate of what the sample VAF would be if the host were exactly
#'   heterozygous
#' @md
#' @examples
#' deconvolute(90, 15, 71, 8, 0.03670243, 0.59, FALSE)
#' @export
deconvolute <- function(tumour_ref, tumour_alt, host_ref, host_alt, logr, purity, plot=FALSE) {
    # All the numerical work is done in find_mle. the rest of this function
    # is for convenient presentation of the result.
    .Deprecated("fast_estimate_tumour_vaf")
    return ()
    result <- find_mle(tumour_ref, tumour_alt, host_ref, host_alt, logr, purity, hess = TRUE)

    covariance <- solve(result$hessian)
    counts <- count_matrix(tumour_ref, tumour_alt, result$par)

    if (plot) {
        # Compute objective over a grid
        X <- seq(0.01, 0.99, length.out = 99)
        Y <- seq(0.01, 0.99, length.out = 99)
        g <- expand.grid(X, Y)
        z <- apply(g, 1, objective)
        g$z <- z

        # Maps objective values to colours: optimum=red
        colours <- colorRamp(c("red","white", "blue"), bias=5, space = "rgb")((z - min(z)) / (max(z - min(z))))
        colours <- apply(colours, 1, function(x) do.call(rgb, as.list(x/255)))

        layout(matrix(c(1,1,1,2), nrow=1))
        plot(g[,1], g[,2], pch = 15, col = colours, asp = TRUE,
             ylab = "P(source=Host|read=Ref)", xlab = "P(source=H|read=Alt)",
             main = "Deconvolution result", cex=2,
             xlim = c(0,1), ylim = c(0,1))
        points(result$par[1], result$par[2], pch = 4, cex = 5, lwd=2)

        # Can use points in `ell` to construct CI around tumour and host vaf estimates
        ell<-ellipse(result$par, abs(covariance), newplot = FALSE, lty=1, lwd=2, npoints=100)
        mats <- apply(ell, 1, count_matrix, ref=tumour_ref, alt=tumour_alt)
        tvafs <- apply(mats, 2, function(x) x[2] / x[3])
        hvafs <- apply(mats, 2, function(x) x[5] / x[6])
        hvafs_min <- which.min(hvafs)
        hvafs_max <- which.max(hvafs)
        tvafs_min <- which.min(tvafs)
        tvafs_max <- which.max(tvafs)

        points(ell[hvafs_min, 1], ell[hvafs_min, 2], pch = 1, col = "goldenrod")
        points(ell[hvafs_max, 1], ell[hvafs_max, 2], pch = 4, col = "goldenrod")
        points(ell[tvafs_min, 1], ell[tvafs_min, 2], pch = 1, col = "dodgerblue")
        points(ell[tvafs_max, 1], ell[tvafs_max, 2], pch = 4, col = "dodgerblue")

        orig_vaf <- tumour_alt / (tumour_alt+tumour_ref)
        xs <- 1:3
        ys <- c(tumour_alt / (tumour_alt+tumour_ref), counts[2,1]/counts[3,1], counts[2,2]/counts[3,2])
        error_bars <- matrix(c(qbeta(0.025, tumour_alt, tumour_ref), qbeta(0.975, tumour_alt, tumour_ref),
                               tvafs[tvafs_min], tvafs[tvafs_max],
                               hvafs[hvafs_min], hvafs[hvafs_max]), nrow = 2, ncol = 3)
        plot(xs, ys, main = "VAF estimates", xlab=NA, ylab = "VAF", xaxt='n', type = "n", xlim = c(0.75, 3.25), ylim=c(0,1))#range(error_bars) * c(1/1.1, 1.1))
        rect(0.975, error_bars[1,1], 1.025, error_bars[2,1], col = "seagreen", border=NA)
        rect(1.975, error_bars[1,2], 2.025, error_bars[2,2], col = "dodgerblue", border = NA)
        rect(2.975, error_bars[1,3], 3.025, error_bars[2,3], col = "goldenrod", border = NA)
        abline(h=ys, col = c("seagreen", "dodgerblue", "goldenrod"), lty=3, lwd=0.5)
        points(xs, ys, pch = 21, cex = 2, bg = c("seagreen", "dodgerblue", "goldenrod"), lwd = 2)
        axis(1, at=1:3, labels = c("Original", "Pure tumour", "Pure host"))
    }

    list(par=result$par, cov=covariance, counts=counts,
         pure_tumour_vaf=counts[2,1] / counts[3,1],
         pure_host_vaf=counts[2,2] / counts[3,2],
         sample_vaf_with_het_host=(counts[2,1] + counts[3,2] / 2) / counts[3,3])
}


### Fast alternatives

#' Fast deconvolution of read count contingency table
#' @param total_readdepth Tumour sample total read depth
#' @param alt_readdepth Tumour sample alt read depth
#' @param logr Tumour sample logR
#' @param host_total_readdepth Host sample total read depth
#' @param host_alt_readdepth Host sample alt read depth
#' @param purity Estimated purity of tumour sample
.estimate_contingency_table <- function(total_readdepth,
                                        alt_readdepth,
                                        logr,
                                        host_total_readdepth,
                                        host_alt_readdepth,
                                        purity,
                                        ploidy,
                                        host_ploidy = 2) {
    # Aim: find values for the empty cells (a), (b), (c), (d), (e), (f)
    # in the read count contingency table below. A and T are observed values, and R = T - A.
    # The other values need to be estimated, subject to the constraint that they are all
    # non-negative.
    # (c) (referred to throughout as variable K) can be estimated as T*p_1, where p_1 is the
    # probability that a read came from host material in the mixed tumour-host source. p_1
    # is estimated using the tumour sample's logR (logr), and the host sample's VAF (hvaf).
    # (b) (variable L) is estimated from (c) as K*p_2, where p_2 is the probability that
    # a host-derived read carries the Alt allele, not the Ref allele. p_2 is equal to hvaf.
    # As the table has 2 degrees of freedom, once values are estimated for (b) and (c), the
    # other missing values are filled in trivially.

    # Contingency table :
    #'         |  Ref  |  Alt  |
    #'  -------|-------|-------|-------
    #'   Host  |  (a)  |  (b)  |  (c)
    #'  Tumour |  (d)  |  (e)  |  (f)
    #'  -------|-------|-------|-------
    #'         |   R   |   A   |   T

    hvaf[is.nan(hvaf)] <- 0
    r <- 2^logr
    T <- total_readdepth
    A <- alt_readdepth

    # probability that a read comes from the host (P(R=H))
    p_1 <- 1 - prob_read_came_from_tumour(logr, purity, ploidy, host_ploidy)

    # probability that a host read is an Alt allele (P(R=A|R=H))
    p_2 <- host_alt_readdepth / host_total_readdepth

    # Estimate K (number of host reads) and L (number of host reads that are Alt),
    # subject to the constraints A - L >= 0, and (T - K) - (A - L) >= 0

    K <- T*p_1
    L <- K*p_2

    # If constraints are broken, then find the optimal least squares solution to the
    # constrained equation, using Lagrangian multipliers.
    # Least-squares estimate for K and L:
    # f(K, L) = (K - T*p_1)^2 + (L - K*p_2)
    #
    # Constraints:
    # g1(L, a) => A - L - a^2 = 0
    # g2(K, L, b) => (T-K) - (A-L) - b^2 = 0
    #
    # Augmented Lagrangian equation:
    # F(K, L, a, b, λ1, λ2) = f(K, L) + λ1*g1(L, a) + λ2*g2(K, L, b)
    #
    # Optimal constrained solution obtained when all derivatives of F are zero.

    # Broken constraint 1: Too many Alt reads in the host
    ix <- (L > A)

    K[ix] <- (A[ix]*p_2[ix] + T[ix]*p_1[ix]) / (p_2[ix]^2 + 1)
    L[ix] <- A[ix]

    # Broken constraint 2: VAF out of bounds [ignore if T==0 - these will give NA in the index, which breaks]
    iy <- (T > 0) & ((((A-L) / (T-K)) > 1) | (((A-L) / (T-K)) < 0) | is.nan((A-L) / (T-K)))
    denom <- (p_2[iy]^2 - 2 * p_2[iy] + 2)
    K[iy] <- (A[iy] * (p_2[iy] - 1) + T[iy] * (p_1[iy] - p_2[iy] + 1)) / denom
    L[iy] <- (A[iy] * (p_2[iy]^2 - p_2[iy] + 1) + T[iy] * (p_1[iy] - p_2[iy]^2 + p_2[iy] - 1) ) / denom

    # Contingency table :
    #'     |  Ref  |  Alt  |
    #'  ---|-------|-------|-------
    #'   H |  K-L  |   L   |   K
    #'   T |T-K-A+L|  A-L  |  T-K
    #'  ---|-------|-------|-------
    #'     |   R   |   A   |   T

    return(list(K=K, L=L))
}

#' Probability that a read came from the tumour
#' @param logr Tumour sample logR
#' @param purity Tumour sample purity
#' @param ploidy Tumour ploidy (ploidy of pure tumour)
#' @param host_ploidy Host ploidy estimate (almost always will be 2)
#' @export
#' Derivation:
#' Expected proportion of tumour reads in a mixed sample, where Nt is
#' underlying tumour copy number state, Nh is host copy state, and p
#' is purity:
#'   P(read=Tumour) = p * Nt / (p * Nt + (1 - p) * Nh) (1)
#'
#' Tumour copynumber state is estimated as:
#'   Nt = (R * Nh * (p*ψt + (1-p) * ψh) - ψh * (1 - p) * Nh)
#'        --------------------------------------------------  (2)
#'                            p * ψh
#' Therefore,
#'  p * Nt = (R * Nh * (p*ψt + (1-p) * ψh) - ψh * (1 - p) * Nh)
#'           --------------------------------------------------  (3)
#'                             ψh
#' And,
#'  p * Nt + (1 - p) * Nh
#'         = (R * Nh * (p*ψt + (1-p) * ψh)
#'           -----------------------------   (4)
#'                       ψh
#' Substitute (3) and (4) into (1) to obtain the result.
prob_read_came_from_tumour <- function(logr, purity, ploidy, host_ploidy) {
    stopifnot(purity >= 0 & purity <= 1)
    stopifnot(ploidy > 0)
    stopifnot(host_ploidy > 0)
    denom <- 2^logr * (purity * ploidy + (1 - purity) * host_ploidy)
    p <- (denom - host_ploidy * (1 - purity)) / denom
    pmax(0, pmin(1, p))
}

#' Probability that a read came from the host. See `prob_read_came_from_tumour`
#' @param logr Tumour sample logR
#' @param purity Tumour sample purity
#' @param ploidy Tumour ploidy (ploidy of pure tumour)
#' @param host_copynumber Host copy number estimate (usually 2)
#' @param host_ploidy Host ploidy estimate (almost always will be 2)
#' @export
prob_read_came_from_host <- function(logr, purity, ploidy, host_ploidy) {
    1 - prob_read_came_from_tumour(logr, purity, ploidy, host_ploidy)
}

#' Fast approximate estimate of pure_tumour_vaf from a mixed sample
#' @param total_readdepth Tumour sample total read depth
#' @param alt_readdepth Tumour sample alt read depth
#' @param logr Tumour sample logR
#' @param host_total_readdepth Host sample total read depth
#' @param host_alt_readdepth Host sample alt read depth
#' @param purity Estimated purity of tumour sample
#' @param ploidy estimated ploidy of tumour sample
#' @param host_ploidy estimated ploidy of host sample (default = 2)
#' @export
fast_estimate_tumour_vaf <- function(total_readdepth, alt_readdepth, logr, host_total_readdepth, host_alt_readdepth, purity, ploidy, host_ploidy = 2) {
    result <- .estimate_contingency_table(total_readdepth, alt_readdepth, logr, host_total_readdepth, host_alt_readdepth, purity, ploidy, host_ploidy = 2)
    alt_reads <- alt_readdepth - result$L
    total_reads <- total_readdepth - result$K

    eps <- sqrt(.Machine$double.eps)
    alt_reads[abs(alt_reads) < eps] <- 0
    total_reads[abs(total_reads) < eps] <- 0

    vaf <- alt_reads / total_reads
    vaf[is.na(vaf)] <- 0
    pmax(0, pmin(1, vaf))
}

#' Fast approximate estimate of log-odds ratio from a mixed sample,
#' with the host converted to heterozygous diploid
#' @param total_readdepth Tumour sample total read depth
#' @param alt_readdepth Tumour sample alt read depth
#' @param logr Tumour sample logR
#' @param host_total_readdepth Host sample total read depth
#' @param host_alt_readdepth Host sample alt read depth
#' @param purity Estimated purity of tumour sample
#' @param ploidy Estimated ploidy of the tumour sample
#' @param host_ploidy Estimated ploidy of the host sample (default = 2)
#' @param pseudocount Avoid infinite log-odds by adding a pseudocount to each
#'     cell of the contingency table. Default = 0.5 (same as Facets).
#' @export
fast_estimate_tumour_logodds <- function(total_readdepth, alt_readdepth, logr, host_total_readdepth, host_alt_readdepth, purity, ploidy, host_ploidy = 2, pseudocount=0.5) {
    result <- .estimate_contingency_table(total_readdepth, alt_readdepth, logr, host_total_readdepth, host_alt_readdepth, purity, ploidy, host_ploidy = 2)
    T <- total_readdepth
    A <- alt_readdepth
    L <- result$L
    K <- result$K
    a = K/2 + pseudocount
    b = a
    c = T - K - A + L + pseudocount
    d = A - L + pseudocount
    return (log((a*d) / (b*c)))
}

#' Fast approximate estimate of log-odds ratio, and its variance, from a mixed sample,
#' with the host converted to heterozygous diploid.
#' @param total_readdepth Tumour sample total read depth
#' @param alt_readdepth Tumour sample alt read depth
#' @param logr Tumour sample logR
#' @param host_total_readdepth Host sample total read depth
#' @param host_alt_readdepth Host sample alt read depth
#' @param purity Estimated purity of tumour sample
#' @param ploidy Estimated ploidy of the tumour sample
#' @param host_ploidy Estimated ploidy of the host sample (default = 2)
#' @param pseudocount Avoid infinite log-odds by adding a pseudocount to each
#'     cell of the contingency table. Default = 0.5 (same as Facets).
#' @export
fast_estimate_tumour_logodds_and_variance <- function(total_readdepth, alt_readdepth, logr, host_total_readdepth, host_alt_readdepth, purity, ploidy, host_ploidy = 2, pseudocount=0.5) {
    result <- .estimate_contingency_table(total_readdepth, alt_readdepth, logr, host_total_readdepth, host_alt_readdepth, purity, ploidy, host_ploidy = 2)
    T <- total_readdepth
    A <- alt_readdepth
    L <- result$L
    K <- result$K
    a = K/2 + pseudocount
    b = a
    c = T - K - A + L + pseudocount
    d = A - L + pseudocount
    logodds <- (log((a*d) / (b*c)))
    var_logodds <- 1/a + 1/b + 1/c + 1/d
    return (list(logodds = logodds, var_logodds = var_logodds))
}

#' Fast approximate estimate of pure_tumour_vaf from a mixed sample
#' @param total_readdepth Tumour sample total read depth
#' @param alt_readdepth Tumour sample alt read depth
#' @param logr Tumour sample logR
#' @param host_total_readdepth Host sample total read depth
#' @param host_alt_readdepth Host sample alt read depth
#' @param purity Estimated purity of tumour sample
#' @param ploidy Estimated ploidy of the tumour sample
#' @param host_ploidy Estimated ploidy of the host sample (default = 2)
#' @export
fast_estimate_contingency_table <- function(total_readdepth, alt_readdepth, logr, host_total_readdepth, host_alt_readdepth, purity, ploidy, host_ploidy = 2) {
    result <- .estimate_contingency_table(total_readdepth, alt_readdepth, logr, host_total_readdepth, host_alt_readdepth, purity, ploidy, host_ploidy)
    T <- total_readdepth
    A <- alt_readdepth
    L <- result$L
    K <- result$K
    matrix(c(K - L, T - K - A + L, T - A, L, A - L, A, K, T - K, T), nrow = 3)
}

.compound_vaf_correct <- function(total_readddepth,
                                  alt_readddepth,
                                  logr,
                                  host_total_readdepth,
                                  host_alt_readdepth,
                                  purity,
                                  ploidy,
                                  host_ploidy = 2,
                                  use_pseudocount = TRUE) {
    T <- total_readddepth
    A <- alt_readddepth
    p_1 <- prob_read_came_from_host(logr, purity, ploidy, host_ploidy)
    alpha <- host_alt_readdepth
    beta <- host_total_readdepth - host_alt_readdepth

    if (use_pseudocount) {
        alpha <- alpha + 0.5
        beta <- beta + 0.5
    }

    g <- expand.grid(K = 0:T, L = 0:A)
    g <- as.data.table(g)[(L >= A - T + K) & (L <= A)]
    g[, p := dbetabinom.ab(L, K, alpha, beta, log = FALSE) * dbinom(K, T, p_1), by = K]
    g[, p := p / sum(p)]

    expectedL <- g[, sum(L*p)]
    expectedK <- g[, sum(K*p)]
    vaf <- (A - expectedL) / (T - expectedK)
    ifelse(is.na(vaf), 0, vaf)
}

#' @export
compound_vaf_correct <- Vectorize(.compound_vaf_correct,
                                  vectorize.args = c(
                                      "total_readddepth",
                                      "alt_readddepth",
                                      "logr",
                                      "host_total_readdepth",
                                      "host_alt_readdepth"))

