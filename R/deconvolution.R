#' Fast deconvolution of read count contingency table. Used as basis of
#' the exported fast_estimate_* functions
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
    # (c) (referred to throughout as variable K) can be estimated as T*p_host, where p_host is the
    # probability that a read came from host material in the mixed tumour-host source. p_host
    # is estimated using the tumour sample's logR (logr), and the host sample's VAF (hvaf).
    # (b) (variable L) is estimated from (c) as K*p_alt_given_host, where p_alt_given_host is the probability that
    # a host-derived read carries the Alt allele, not the Ref allele. p_alt_given_host is equal to hvaf.
    # As the table has 2 degrees of freedom, once values are estimated for (b) and (c), the
    # other missing values are filled in trivially.

    # Contingency table :
    #'         |  Ref  |  Alt  |
    #'  -------|-------|-------|-------
    #'   Host  |  (a)  |  (b)  |  (c)
    #'  Tumour |  (d)  |  (e)  |  (f)
    #'  -------|-------|-------|-------
    #'         |   R   |   A   |   T

    hvaf <- host_alt_readdepth / host_total_readdepth
    hvaf[is.nan(hvaf)] <- 0
    r <- 2^logr
    T <- total_readdepth
    A <- alt_readdepth

    # probability that a read comes from the host (P(R=H))
    p_host <- prob_read_came_from_host(logr, purity, ploidy, host_ploidy)

    # probability that a host read is an Alt allele (P(R=A|R=H))
    p_alt_given_host <- hvaf

    # Estimate K (number of host reads) and L (number of host reads that are Alt),
    # subject to the constraints A - L >= 0, and (T - K) - (A - L) >= 0

    K <- T*p_host
    L <- K*p_alt_given_host

    # If constraints are broken, then find the optimal least squares solution to the
    # constrained equation, using Lagrangian multipliers.
    # Least-squares estimate for K and L:
    # f(K, L) = (K - T*p_host)^2 + (L - K*p_alt_given_host)
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

    K[ix] <- (A[ix]*p_alt_given_host[ix] + T[ix]*p_host[ix]) / (p_alt_given_host[ix]^2 + 1)
    L[ix] <- A[ix]

    # Broken constraint 2: VAF out of bounds [ignore if T==0 - these will give NA in the index, which breaks]
    iy <- (T > 0) & ((((A-L) / (T-K)) > 1) | (((A-L) / (T-K)) < 0) | is.nan((A-L) / (T-K)))
    denom <- (p_alt_given_host[iy]^2 - 2 * p_alt_given_host[iy] + 2)
    K[iy] <- (A[iy] * (p_alt_given_host[iy] - 1) + T[iy] * (p_host[iy] - p_alt_given_host[iy] + 1)) / denom
    L[iy] <- (A[iy] * (p_alt_given_host[iy]^2 - p_alt_given_host[iy] + 1) + T[iy] * (p_host[iy] - p_alt_given_host[iy]^2 + p_alt_given_host[iy] - 1) ) / denom

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
#' @export
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
    matrix(c(A - L, T - A - K + L, T - K,
             L, K - L, K,
             A, T - A, T),
           nrow = 3, ncol = 3,
           dimnames = list(c("Alt", "Ref", "Total"),
                           c("Tumour", "Host", "Sample")))
}

#' Probability mass function for the beta-binomial distribution
#' @param x Integer: Number of successes
#' @param n Integer: Number of trials
#' @param a Numeric: alpha parameter of the beta distribution
#' @param b Numeric: beta parameter of the beta distribution
#' @returns Numeric: probability of x given n, a, b
#' @export
dbetabinom <- function(x, n, a, b) {
    e1 <- gamma(n + 1) / (gamma(x + 1) * gamma(n - x + 1))
    e2 <- beta(x + a, n - x + b) / beta(a, b)
    e1 * e2
}

#' Calculates the expectation of the tumour VAF, marginalized over
#' all valid contingency tables, conditioned on observed read totals
#' Not an exported function - instead expectation_tumour_vaf is the vectorized form of this
#' @param total_readdepth Integer: Observed total read depth in the tumour sample
#' @param alt_readdepth Integer: Observed number of reads carrying the alt allele in the tumour sample
#' @param logr Numeric: Estimated log2 of the ratio of tumour read depth/host read depth
#' @param host_total_readdepth Integer: Observed total read depth in the matched host sample
#' @param host_alt_readdepth Integer: Observed number of reads carrying the alt allele in the matched host sample
#' @param purity Numeric [0, 1]: Estimated aberrant cell fraction of the tumour sample
#' @param ploidy Numeric: Estimated ploidy of the average tumour cell
#' @param host_ploidy Numeric: Estimated ploidy of a normal cell (default=2)
#' @param use_pseudocount Logical: If TRUE, adds 0.5 to the alpha and beta parameters of the beta-binomial distribution
#' @importFrom data.table as.data.table
.expectation_tumour_vaf <- function(total_readdepth,
                                    alt_readdepth,
                                    logr,
                                    host_total_readdepth,
                                    host_alt_readdepth,
                                    purity,
                                    ploidy,
                                    host_ploidy = 2,
                                    use_pseudocount = TRUE) {

    T <- total_readdepth
    A <- alt_readdepth
    alpha <- host_alt_readdepth
    beta <- host_total_readdepth - alpha

    if (use_pseudocount) {
        alpha <- alpha + 0.5
        beta <- beta + 0.5
    }

    pH <- prob_read_came_from_host(logr, purity, ploidy, 2)

    # Enumerate all possible contingency tables (two degrees of freedom)
    dt <- as.data.table(expand.grid(Ta = 0:A, Tt = 0:T))[Ta <= Tt]
    dt[, L := A - Ta]
    dt[, K := T - Tt]
    # Restrict to valid tables (all elements non-negative integers; row-marginal sums are respected)
    dt <- dt[(L >= A - T + K) & (L >= 0) & (L <= A) & (L <= K)]
    # Compute probability of observing L host alt reads out of K host reads (beta-binomial),
    # weighted by probability of observing K host reads out of T total reads (binomial)

    .probfun <- function(l, k, t, a, b, p) {
        # betabinom
        lgk <- lgamma(k + 1)
        coeff_bb <- lgk - lgamma(l + 1) - lgamma(k - l + 1)
        logp_bb <- lbeta(l + a, k - l + b) - lbeta(a, b)

        # binom
        logp_b <- dbinom(k, t, p, log = TRUE)
        exp(coeff_bb + logp_bb + coeff_b + logp_b)
    }

    dt[, p := .probfun(L, K, T, alpha, beta, pH)]
    # Condition probability on table being valid
    dt[, p := p / sum(p)]

    # Compute VAF for each entry, and return its expectation
    dt[, vaf := ifelse(Tt == 0, 0, Ta / Tt)]
    dt[, sum(p * vaf)]
}

#' Calculates the expectation of the tumour VAF, marginalized over
#' all valid contingency tables, conditioned on observed read totals.
#' Much slower than `fast_estimate_tumour_vaf`
#' @param total_readdepth Integer: Observed total read depth in the tumour sample
#' @param alt_readdepth Integer: Observed number of reads carrying the alt allele in the tumour sample
#' @param logr Numeric: Estimated log2 of the ratio of tumour read depth/host read depth
#' @param host_total_readdepth Integer: Observed total read depth in the matched host sample
#' @param host_alt_readdepth Integer: Observed number of reads carrying the alt allele in the matched host sample
#' @param purity Numeric [0, 1]: Estimated aberrant cell fraction of the tumour sample
#' @param ploidy Numeric: Estimated ploidy of the average tumour cell
#' @param host_ploidy Numeric: Estimated ploidy of a normal cell (default=2)
#' @param use_pseudocount Logical: If TRUE, adds 0.5 to the alpha and beta parameters of the beta-binomial distribution
#' @export
expectation_tumour_vaf <- Vectorize(.expectation_tumour_vaf,
                                  vectorize.args = c(
                                      "total_readdepth",
                                      "alt_readdepth",
                                      "logr",
                                      "host_total_readdepth",
                                      "host_alt_readdepth"))
