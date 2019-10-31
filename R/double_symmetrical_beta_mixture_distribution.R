# Double symmetric Beta
# Functions describing a mixture of two mixtures of two Beta distributions with symmetrical parameters,
# i.e. w3w1 * Beta(A, B) + w3(1-w1) * Beta(B, A) + (1-w3)w2 * Beta(C, D) + (1-w3)(1-w2) * Beta(D, C)
# w1, w2 and w3 are mixture weights:
# w1 and w2 control the inner weights of each symmetrical Beta mixture, and w3 controls the weight
# of the mixture of mixtures

# Why? Allele-specific copy number estimation requires VAF estimated from (pre-copynumber-change-) heterozygous SNVs,
# but data consists of a mixture of homozygous and heterozygous SNVs. A symmetrical Beta is a good model for
# estimating VAF

#' Density function for symmetrical beta mixture
#' @export
ddoublesbeta <- function(x, a1, b1, a2, b2, w1, w2, w3, log = FALSE) {
    if (log) {
        d <- data.frame(e1=log(w1) + log(w3) + dbeta(x, a1, b1, log = TRUE),
                        e2=log(1-w1) + log(w3) + dbeta(x, b1, a1, log = TRUE),
                        e3=log(w2) + log(1-w3) + dbeta(x, a2, b2, log = TRUE),
                        e4=log(1-w2) + log(1-w3) + dbeta(x, b2, a2, log = TRUE))
        # Mixture of two densities is done elementwise, using logSumExp
        return(apply(d, 1, function(r) {
            m <- max(r)
            log(sum(exp(r - m))) + m
        }))
    }
    # Mixture of two densities
    w1 * w3 * dbeta(x, a1, b1) + (1 - w1) * w3 * dbeta(x, b1, a1) +
        w2 * (1 - w3) * dbeta(x, a2, b2) + (1 - w2) * (1 - w3) * dbeta(x, b2, a2)
}

#' PDF function for symmetrical beta mixture
#' @export
pdoublesbeta <- function(q, a1, b1, a2, b2, w1, w2, w3) {
    w1 * w3 * pbeta(q, a1, b1) + (1 - w1) * w3 * pbeta(q, b1, a1) +
        w2 * (1 - w3) * pbeta(q, a2, b2) + (1 - w2) * (1 - w3) * pbeta(q, b2, a2)
}

#' Quantile (CDF) function for symmetrical beta mixture
#' @export
qdoublesbeta <- function(p, a1, b1, a2, b2, w1, w2, w3) {
    G <- function(x) pdoublesbeta(x, a1, b1, a2, b2, w1, w2, w3) - p
    return (uniroot(G, c(0, 1))$root)
}

#' Random variates function for symmetrical beta mixture
#' @export
rdoublesbeta <- function(n, a1, b1, a2, b2, w1, w2, w3) {
    sample(c(rbeta(n*w1*w3, a1, b1),
             rbeta(n*(1-w1)*w3, b1, a1),
             rbeta(n*w2*(1-w3),a2, b2),
             rbeta(n*(1-w2)*(1-w3), b2, a2)))
}

# Derivatives of the pdf
ddoublesbeta_deriv_wrt_a1 <- function(x, a1, b1, a2, b2, w1, w2, w3) {
    -w3*(w1*x^a1*(-x + 1)^b1*(digamma(a1 + 1) - digamma(a1 + b1 + 2)) - w1*x^a1*(-x + 1)^b1*log(x) - x^b1*(w1 - 1)*(-x + 1)^a1*(digamma(a1 + 1) - digamma(a1 + b1 + 2)) + x^b1*(w1 - 1)*(-x + 1)^a1*log(-x + 1))/beta(a1 + 1, b1 + 1)
}

ddoublesbeta_deriv_wrt_a2 <- function(x, a1, b1, a2, b2, w1, w2, w3) {
    (w3 - 1)*(w2*x^a2*(-x + 1)^b2*(digamma(a2 + 1) - digamma(a2 + b2 + 2)) - w2*x^a2*(-x + 1)^b2*log(x) - x^b2*(w2 - 1)*(-x + 1)^a2*(digamma(a2 + 1) - digamma(a2 + b2 + 2)) + x^b2*(w2 - 1)*(-x + 1)^a2*log(-x + 1))/beta(a2 + 1, b2 + 1)
}

ddoublesbeta_deriv_wrt_b1 <- function(x, a1, b1, a2, b2, w1, w2, w3) {
    -w3*(w1*x^a1*(-x + 1)^b1*(digamma(b1 + 1) - digamma(a1 + b1 + 2)) - w1*x^a1*(-x + 1)^b1*log(-x + 1) - x^b1*(w1 - 1)*(-x + 1)^a1*(digamma(b1 + 1) - digamma(a1 + b1 + 2)) + x^b1*(w1 - 1)*(-x + 1)^a1*log(x))/beta(a1 + 1, b1 + 1)
}

ddoublesbeta_deriv_wrt_b2 <- function(x, a1, b1, a2, b2, w1, w2, w3) {
    (w3 - 1)*(w2*x^a2*(-x + 1)^b2*(digamma(b2 + 1) - digamma(a2 + b2 + 2)) - w2*x^a2*(-x + 1)^b2*log(-x + 1) - x^b2*(w2 - 1)*(-x + 1)^a2*(digamma(b2 + 1) - digamma(a2 + b2 + 2)) + x^b2*(w2 - 1)*(-x + 1)^a2*log(x))/beta(a2 + 1, b2 + 1)
}

ddoublesbeta_deriv_wrt_w1 <- function(x, a1, b1, a2, b2, w1, w2, w3) {
    w3*(x^a1*(-x + 1)^b1 - x^b1*(-x + 1)^a1)/beta(a1 + 1, b1 + 1)
}

ddoublesbeta_deriv_wrt_w2 <- function(x, a1, b1, a2, b2, w1, w2, w3) {
    -(w3 - 1)*(x^a2*(-x + 1)^b2 - x^b2*(-x + 1)^a2)/beta(a2 + 1, b2 + 1)
}

ddoublesbeta_deriv_wrt_w3 <- function(x, a1, b1, a2, b2, w1, w2, w3) {
    ((w1*x^a1*(-x + 1)^b1 + x^b1*(-w1 + 1)*(-x + 1)^a1)*beta(a2 + 1, b2 + 1) + (-w2*x^a2*(-x + 1)^b2 + x^b2*(w2 - 1)*(-x + 1)^a2)*beta(a1 + 1, b1 + 1))/(beta(a1 + 1, b1 + 1)*beta(a2 + 1, b2 + 1))
}


#' Log-likelihood of a sample conditioned on these parameters, assuming symmetrical beta mixture model
#' @export
lldoublesbeta <- function(xs, a1, b1, a2, b2, w1, w2, w3) {
    sum(ddoublesbeta(xs, a1, b1, a2, b2, w1, w2, w3, log = TRUE))
}

# Derivatives of the likelihood
lldoublesbeta_deriv_wrt_a1 <- function(xs, a1, b1, a2, b2, w1, w2, w3) {
    sum(ddoublesbeta_deriv_wrt_a1(xs, a1, b1, a2, b2, w1, w2, w3) / ddoublesbeta(xs, a1, b1, a2, b2, w1, w2, w3))
}

lldoublesbeta_deriv_wrt_a2 <- function(xs, a1, b1, a2, b2, w1, w2, w3) {
    sum(ddoublesbeta_deriv_wrt_a2(xs, a1, b1, a2, b2, w1, w2, w3) / ddoublesbeta(xs, a1, b1, a2, b2, w1, w2, w3))
}

lldoublesbeta_deriv_wrt_b1 <- function(xs, a1, b1, a2, b2, w1, w2, w3) {
    sum(ddoublesbeta_deriv_wrt_b1(xs, a1, b1, a2, b2, w1, w2, w3) / ddoublesbeta(xs, a1, b1, a2, b2, w1, w2, w3))
}

lldoublesbeta_deriv_wrt_b2 <- function(xs, a1, b1, a2, b2, w1, w2, w3) {
    sum(ddoublesbeta_deriv_wrt_b2(xs, a1, b1, a2, b2, w1, w2, w3) / ddoublesbeta(xs, a1, b1, a2, b2, w1, w2, w3))
}

lldoublesbeta_deriv_wrt_w3 <- function(xs, a1, b1, a2, b2, w1, w2, w3) {
    sum(ddoublesbeta_deriv_wrt_w3(xs, a1, b1, a2, b2, w1, w2, w3) / ddoublesbeta(xs, a1, b1, a2, b2, w1, w2, w3))
}

lldoublesbeta_deriv_wrt_w1 <- function(xs, a1, b1, a2, b2, w1, w2, w3) {
    sum(ddoublesbeta_deriv_wrt_w1(xs, a1, b1, a2, b2, w1, w2, w3) / ddoublesbeta(xs, a1, b1, a2, b2, w1, w2, w3))
}

lldoublesbeta_deriv_wrt_w2 <- function(xs, a1, b1, a2, b2, w1, w2, w3) {
    sum(ddoublesbeta_deriv_wrt_w2(xs, a1, b1, a2, b2, w1, w2, w3) / ddoublesbeta(xs, a1, b1, a2, b2, w1, w2, w3))
}


finite_diff <- function(fn, params, param_index, h = 0.0001) {
    params_1 <- as.list(params)
    params_2 <- as.list(params)
    for (i in param_index) {
        params_1[[i]] <- params_1[[i]] - h
        params_2[[i]] <- params_2[[i]] + h
    }
    (do.call(fn, params_2) - do.call(fn, params_1)) / (2*h)
}

#' Estimate parameters of symmetrical Beta distribution by maximum likelihood
#' (assumes mixture weight fixed to 0.5)
#' @importFrom "nloptr" nl.opts
#' @export
est_doublesbeta <- function(data, startval = c(1, 10, 1, 10, 0.5, 0.5, 0.5)) {
    objective <- function(param) {
        -lldoublesbeta(data, param[1], param[2], param[3], param[4], param[5], param[6], param[7])
    }

    gradient <- function(param) {
        -c(
            lldoublesbeta_deriv_wrt_a1(data, param[1], param[2], param[3], param[4], param[5], param[6], param[7]),
            lldoublesbeta_deriv_wrt_b1(data, param[1], param[2], param[3], param[4], param[5], param[6], param[7]),
            lldoublesbeta_deriv_wrt_a2(data, param[1], param[2], param[3], param[4], param[5], param[6], param[7]),
            lldoublesbeta_deriv_wrt_b2(data, param[1], param[2], param[3], param[4], param[5], param[6], param[7]),
            lldoublesbeta_deriv_wrt_w1(data, param[1], param[2], param[3], param[4], param[5], param[6], param[7]),
            lldoublesbeta_deriv_wrt_w2(data, param[1], param[2], param[3], param[4], param[5], param[6], param[7]),
            lldoublesbeta_deriv_wrt_w3(data, param[1], param[2], param[3], param[4], param[5], param[6], param[7])
        )
    }

    # optim(startval, objective, gradient, method = "L-BFGS-B", control = list(trace=100),
    #       lower = c(0.01,0.01,0.01,0.01,0.01,0.01,0.01),
    #       upper = c(999,999,999,999,0.99,0.99,0.99))

    nloptr::auglag(startval, objective, gradient,
                  lower = c(0.01,0.01,0.01,0.01,0.01,0.01,0.01),
                  upper = c(999,999,999,999,0.99,0.99,0.99),
                  hin = function(param) (param[2] - param[1]) * (param[4] - param[3]), control = list(xtol_rel = 1e-8, maxeval = 2000))
}


#' Estimate parameters of symmetrical Beta distribution by maximum likelihood
#' (assumes mixture weight fixed to 0.5)
#' @export
est_doublesbeta2 <- function(data, startval = c(1, 10, 1, 10, 0.5)) {
    objective <- function(param) {
        -lldoublesbeta(data, param[1], param[2], param[3], param[4], 0.5, 0.5, param[5])
    }

    gradient <- function(param) {
        -c(
            lldoublesbeta_deriv_wrt_a1(data, param[1], param[2], param[3], param[4], 0.5, 0.5, param[5]),
            lldoublesbeta_deriv_wrt_b1(data, param[1], param[2], param[3], param[4], 0.5, 0.5, param[5]),
            lldoublesbeta_deriv_wrt_a2(data, param[1], param[2], param[3], param[4], 0.5, 0.5, param[5]),
            lldoublesbeta_deriv_wrt_b2(data, param[1], param[2], param[3], param[4], 0.5, 0.5, param[5]),
            lldoublesbeta_deriv_wrt_w3(data, param[1], param[2], param[3], param[4], 0.5, 0.5, param[5])
        )
    }

    # optim(startval, objective, gradient, method = "L-BFGS-B",
    #       lower = c(0.01,0.01,0.01,0.01,0.01,0.01,0.01),
    #       upper = c(999,999,999,999,0.99,0.99,0.99))

    nloptr::slsqp(startval, objective, gradient,
                  lower = c(0.01,0.01,0.01,0.01, 0.01),
                  upper = c(999,999,999,999,0.99),
                  hin = function(param) (param[2] - param[1]) * (param[4] - param[3]),control = nl.opts(list(xtol_rel = 1e-8, maxeval = 2000)))
}

# (Debug function)
.add_tangent_line <- function(a1, b1, a2, b2, w1, w2, w3, data) {
    y <- lldoublesbeta(data, a1, b1, a2, b2, w1, w2, w3)
    grad <- lldoublesbeta_deriv_wrt_a1(data, a1, b1, a2, b2, w1, w2, w3)
    i <- y - grad*a1
    print(grad)
    abline(i, grad)
}


