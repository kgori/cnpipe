#' Log-likelihood for deconvoluting host and tumour read counts
#' @param PHa = parameter probability read is Host given read is A (P(H|A))
#' @param PHb = parameter probability read is Host given read is B (P(H|B))
#' @param n = normal read count in tumour
#' @param v = variant read count in tumour
#' @param logr = log ratio of tumour depth to host depth
#' @param pur = tumour purity
#' @param a = Expected number of host reads
#' @param b = Expected number of tumour reads
loglik <- function(PHa, PHb, n, v, logr, pur, a, b) {
    PHan <- PHa * n
    PHbv <- PHb * v
    PHan_add_PHbv <- PHan + PHbv
    n_add_v <- n + v
    r <- exp(logr)

    # Expected number of tumour reads
    Et <- ((n + v - 1) * r * pur) / (r * pur + 1 - pur)

    # Expected number of host reads
    Eh <- (1-pur) * (n + v - 1) / (r * pur + 1 - pur)

    (Eh - 1) * log((PHan_add_PHbv) / n_add_v) + (Et - 1) * log((n_add_v - PHan_add_PHbv)/ n_add_v) +
        (a - 1) * log(PHan / (PHan_add_PHbv)) + (b - 1) * log(PHbv / (PHan_add_PHbv)) - lbeta(Et, Eh) - lbeta(a, b)
}

#' Gradient of log-likelihood for deconvoluting host and tumour read counts w.r.t PHa and PHb
#' @param PHa = parameter probability read is Host given read is A (P(H|A))
#' @param PHb = parameter probability read is Host given read is B (P(H|B))
#' @param n = normal read count in tumour
#' @param v = variant read count in tumour
#' @param logr = log ratio of tumour depth to host depth
#' @param pur = tumour purity
#' @param a = Expected number of host reads
#' @param b = Expected number of tumour reads
dloglik <- function(PHa, PHb, n, v, logr, pur, a, b) {
    r <- exp(logr)

    # Expected total number of tumour reads
    ETt <- ((n + v - 1) * r * pur) / (r * pur + 1 - pur)

    # Expected number of host reads
    EHt <- (1-pur) * (n + v - 1) / (r * pur + 1 - pur)

    # Expected number of host A reads
    EHa <- PHa*n

    # Expected number of host B reads
    EHb <- PHb*v

    # Define terms used repeatedly to save multiplications
    EHa_sq <- EHa*EHa
    EHb_sq <- EHb*EHb
    EHa_EHb <- EHa*EHb
    denom_term <- (EHa_sq + 2*EHa_EHb - EHa*n - EHa*v + EHb_sq - EHb*n - EHb*v)

    dPHa <- (EHt*EHa_EHb + EHt*EHa_sq - EHt*EHa*v - EHt*EHa*n + ETt*EHa_EHb + ETt*EHa_sq - EHa_sq*b - EHa_sq + EHa_EHb*a - EHa_EHb*b - 2*EHa_EHb + EHa*b*n + EHa*b*v + EHb_sq*a - EHb_sq - EHb*a*n - EHb*a*v + EHb*n + EHb*v)/(PHa*denom_term)
    dPHb <- (EHt*EHa_EHb + EHt*EHb_sq - EHt*EHb*n - EHt*EHb*v + ETt*EHa_EHb + ETt*EHb_sq + EHa_sq*b - EHa_sq - EHa_EHb*a + EHa_EHb*b - 2*EHa_EHb - EHa*b*n - EHa*b*v - EHb_sq*a - EHb_sq + EHb*a*n + EHb*a*v + EHa*v + EHa*n)/(PHb*denom_term)
    c(dPHa, dPHb)
}

#' Construct a matrix of read counts given the marginal REF and ALT totals,
#' and the conditional probabilities that a read came from the host, given its
#' status as an A (ref) or a B (alt) read
#' Table:
#'
#'     |   T   |   H    | total
#'  ---|-------|--------|------
#'   A |   .   | p(H|A) |  REF
#'   B |   .   | p(H|B) |  ALT
#'  ---|-------|--------|------
#'     |Σ(.|T) | Σ(.|H) |REF+ALT
#'
#' @param ref = Integer count of REF reads observed in (contaminated) tumour sample
#' @param alt = Integer count of ALT reads observed in tumour sample
#' @param p = Vector of conditional probabilities [p(H|A), p(H|B)]
count_matrix <- function(ref, alt, p) {
    m <- matrix(0, nrow = 3, ncol = 3)
    m[,3] <- c(ref, alt, ref+alt)
    m[1,2] <- p[1] * ref
    m[2,2] <- p[2] * alt
    m[1,1] <- ref - m[1,2]
    m[2,1] <- alt - m[2,2]
    m[3,1] <- sum(m[1:2,1])
    m[3,2] <- sum(m[1:2,2])
    m
}

#' Estimate the host and tumour components from a host-contaminated tumour sample
#' @param tumour_ref - number of REF reads observed in host-contaminated tumour, at current site
#' @param tumour_alt - number of ALT reads observed in host-contaminated tumour, at current site
#' @param host_ref - number of REF reads observed in matched host, at current site
#' @param host_alt - number of ALT reads observed in matched host, at current site
#' @param logr - log of ratio between tumour coverage and matched host coverage, at current site
#' @param purity - Estimate of tumour purity (aka aberrant cell count)
#' @param plot - Should a plot be drawn?
#' @importFrom "mixtools" ellipse
#' @return list:
#' par = optimised parameters (p(H|A), p(H|B))
#' cov = estimated covariance of optimised parameters
#'       (Fisher information - inverse of Hessian)
#' counts = Matrix of inferred read counts, divided by Host/Tumour, A/B-allele
#' pure_tumour_vaf = Estimate of tumour VAF in absence of host
#' pure_host_vaf = Estimate of host VAF in absence of tumour
#' @export
estimate <- function(tumour_ref, tumour_alt, host_ref, host_alt, logr, purity, plot=FALSE) {
    pseudocount <- 0.5
    n <- tumour_ref + pseudocount
    v <- tumour_alt + pseudocount
    a <- host_ref + pseudocount
    b <- host_alt + pseudocount

    objective <- function(x) {
        -loglik(x[1], x[2], n, v, logr, purity, a, b)
    }

    gradient <- function(x) {
        -dloglik(x[1], x[2], n, v, logr, purity, a, b)
    }

    result <- optim(c(0.5, 0.5), objective, gr = gradient, method = "L-BFGS-B",
                    lower = c(1e-6, 1e-6),
                    upper = c(1 - 1e-6, 1 - 1e-6),
                    hessian = TRUE)
    covariance <- solve(result$hessian)
    counts <- count_matrix(tumour_ref, tumour_alt, result$par)

    if (plot) {

        X <- seq(0.01, 0.99, length.out = 99)
        Y <- seq(0.01, 0.99, length.out = 99)
        g <- expand.grid(X, Y)
        z <- apply(g, 1, objective)
        g$z <- z

        colours <- colorRamp(c("red","white", "blue"), bias=5, space = "rgb")((z - min(z)) / (max(z - min(z))))
        colours <- apply(colours, 1, function(x) do.call(rgb, as.list(x/255)))

        layout(matrix(c(1,1,1,2), nrow=1))
        plot(g[,1], g[,2], pch = 15, col = colours, asp = TRUE,
             ylab = "P(source=Host|read=Ref)", xlab = "P(source=H|read=Alt)",
             main = "Deconvolution result", cex=2,
             xlim = c(0,1), ylim = c(0,1))
        points(result$par[1], result$par[2], pch = 4, cex = 5, lwd=2)

        # Can use points in `ell` to construct CI around tumour and host vaf estimates
        ell<-mixtools::ellipse(result$par, abs(covariance), newplot = FALSE, lty=1, lwd=2, npoints=100)
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
         pure_host_vaf=counts[2,2] / counts[3,2])
}

# Uses `estimate` to calculate adjusted tumour vaf as it would be if all host contamination were heterozygous
estimate_tumour_vaf_and_correct_to_het <- function(tumour_ref, tumour_alt, host_ref, host_alt, logr, purity) {
    counts <- estimate(tumour_ref, tumour_alt, host_ref, host_alt, logr, purity, plot = FALSE)$counts
    (counts[2,1] + counts[3,2] / 2) / counts[3,3]
}

# Uses `estimate` to calculate adjusted tumour vaf as it would be if there were no host contamination
estimate_tumour_vaf <- function(tumour_ref, tumour_alt, host_ref, host_alt, logr, purity) {
    estimate(tumour_ref, tumour_alt, host_ref, host_alt, logr, purity, plot = FALSE)$pure_tumour_vaf
}

#' Uses `estimate` to calculate adjusted tumour vaf as it would be if there were no host contamination
#' (vectorised over indicated parameters)
#' @param tumour_ref (vectorised) - number of REF reads observed in host-contaminated tumour, at current site
#' @param tumour_alt (vectorised) - number of ALT reads observed in host-contaminated tumour, at current site
#' @param host_ref (vectorised) - number of REF reads observed in matched host, at current site
#' @param host_alt (vectorised) - number of ALT reads observed in matched host, at current site
#' @param logr (vectorised) - log of ratio between tumour coverage and matched host coverage, at current site
#' @param purity (NOT vectorised) - Estimate of tumour purity (aka aberrant cell count)
#' @return vector of VAF values
#' @export
v.estimate_tumour_vaf <-
    Vectorize(estimate_tumour_vaf,
              vectorize.args = c("tumour_ref", "tumour_alt", "host_ref", "host_alt", "logr"))

#' Uses `estimate` to calculate adjusted tumour vaf as it would be if all host contamination were heterozygous
#' (vectorised over indicated parameters)
#' @param tumour_ref (vectorised) - number of REF reads observed in host-contaminated tumour, at current site
#' @param tumour_alt (vectorised) - number of ALT reads observed in host-contaminated tumour, at current site
#' @param host_ref (vectorised) - number of REF reads observed in matched host, at current site
#' @param host_alt (vectorised) - number of ALT reads observed in matched host, at current site
#' @param logr (vectorised) - log of ratio between tumour coverage and matched host coverage, at current site
#' @param purity (NOT vectorised) - Estimate of tumour purity (aka aberrant cell count)
#' @return vector of VAF values
#' @export
v.estimate_tumour_vaf_and_correct_to_het <-
    Vectorize(estimate_tumour_vaf_and_correct_to_het,
              vectorize.args = c("tumour_ref", "tumour_alt", "host_ref", "host_alt", "logr"))
