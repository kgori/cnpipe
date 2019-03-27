# Functions describing a mixture of two Beta distributions with symmetrical parameters,
# i.e. 0.5 * Beta(A, B) + 0.5 * Beta(B, A)

#' Density function for symmetrical beta mixture
#' @export
dsbeta <- function(x, a, b, w, log = FALSE) {
    if (log) {
        d <- data.frame(e1=log(w) + dbeta(x, a, b, log = TRUE),
                        e2=log(1-w) + dbeta(x, b, a, log = TRUE))
        # Mixture of two densities is done elementwise, using logSumExp
        return(apply(d, 1, function(r) {
            m <- max(r)
            log(sum(exp(r - m))) + m
        }))
    }
    # Mixture of two densities
    w * dbeta(x, a, b) + (1 - w) * dbeta(x, b, a)
}

#' PDF function for symmetrical beta mixture
#' @export
psbeta <- function(q, a, b, w) {
    w * pbeta(q, a, b) + (1 - w) * pbeta(q, b, a)
}

#' Quantile (CDF) function for symmetrical beta mixture
#' @export
qsbeta <- function(p, a, b, w) {
    G <- function(x) psbeta(x, a, b, w) - p
    return (uniroot(G, c(0, 1))$root)
}

#' Random variates function for symmetrical beta mixture
#' @export
rsbeta <- function(n, a, b, w) {
    c(rbeta(n*w, a, b), rbeta(n*(1-w), b, a))
}


# Derivatives of the likelihood
dsbeta_deriv_wrt_x <- function(x, a, b, w) {
    numer1 <- w * x^(a - 1) * (1 - b) * (1 - x)^(b - 1)
    numer2 <- w * x^(a - 1) * (1 - a) * (1 - x)^(b - 1)
    numer3 <- (1 - w) * x^(b - 1) * (1 - a) * (1 - x)^(a - 1)
    numer4 <- (1 - w) * x^(b - 1) * (1 - b) * (1 - x)^(a - 1)
    denom1 <- (1 - x) * beta(a, b)
    denom2 <- x * beta(a, b)
    numer1 / denom1 - numer2 / denom2 + numer3 / denom1 - numer4 / denom2
}

dsbeta_deriv_wrt_a <- function(x, a, b, w) {
    elem1 <- w * x^(a - 1) * (1 - x)^(b - 1) * (log(x) + psigamma(a + b) - psigamma(a))
    elem2 <- (1 - w) * x^(b - 1) * (1 - x)^(a - 1) * (log(1 - x) + psigamma(a + b) - psigamma(a))
    return((elem1 + elem2) / beta(a, b))
}

dsbeta_deriv_wrt_b <- function(x, a, b, w) {
    elem1 <- w * x^(a - 1) * (1 - x)^(b - 1) * (log(1 - x) - psigamma(b) + psigamma(a + b))
    elem2 <- (1 - w) * x^(b - 1) * (1 - x)^(a - 1) * (log(x) - psigamma(b) + psigamma(a + b))
    return((elem1 + elem2) / beta(a, b))
}

dsbeta_deriv_wrt_w <- function(x, a, b, w) {
    elem1 <- x^(a - 1) * (1 - x)^(b - 1)
    elem2 <- x^(b - 1) * (1 - x)^(a - 1)
    return((elem1 - elem2) / beta(a, b))
}


#' Log-likelihood of a sample conditioned on these parameters, assuming symmetrical beta mixture model
#' @export
llsbeta <- function(xs, a, b, w) {
    sum(dsbeta(xs, a, b, w, log = TRUE))
}

llsbeta_deriv_wrt_a <- function(xs, a, b, w) {
    sum(dsbeta_deriv_wrt_a(xs, a, b, w) / dsbeta(xs, a, b, w))
}

llsbeta_deriv_wrt_b <- function(xs, a, b, w) {
    sum(dsbeta_deriv_wrt_b(xs, a, b, w) / dsbeta(xs, a, b, w))
}

llsbeta_deriv_wrt_w <- function(xs, a, b, w) {
    sum(dsbeta_deriv_wrt_w(xs, a, b, w) / dsbeta(xs, a, b, w))
}

dbeta_deriv_wrt_x <- function(x, a, b) {
    betafn_ab <- beta(a, b)
    numer1 <- x^(a - 1) * (1 - b) * (1 - x)^(b - 1)
    denom1 <- (1 - x) * betafn_ab
    numer2 <- x^(a - 1) * (a - 1) * (1 - x)^(b - 1)
    denom2 <- x * betafn_ab
    (numer1 / denom1) + (numer2 / denom2)
}

dbeta_deriv_wrt_a <- function(x, a, b) {
    (x^(a - 1) * (1 - x)^(b - 1) * (log(x) + psigamma(a + b) - psigamma(a))) / beta(a, b)
}

dbeta_deriv_wrt_b <- function(x, a, b) {
    (x^(a - 1) * (1 - x)^(b - 1) * (log(1 - x) + psigamma(a + b) - psigamma(b))) / beta(a, b)
}


finite_diff <- function(fn, params, param_index, h = 0.0001) {
    params_1 <- params
    params_2 <- params
    for (i in param_index) {
        params_1[[i]] <- params_1[[i]] - h
        params_2[[i]] <- params_2[[i]] + h
    }
    (do.call(fn, params_2) - do.call(fn, params_1)) / (2*h)
}

#' Estimate parameters of symmetrical Beta distribution by maximum likelihood
#' (assumes mixture weight fixed to 0.5)
#' @export
est_sbeta <- function(data, startval = c(10, 10)) {
    objective <- function(param) {
        -llsbeta(data, param[1], param[2], 0.5)
    }

    gradient <- function(param) {
        -c(
            llsbeta_deriv_wrt_a(data, param[1], param[2], 0.5),
            llsbeta_deriv_wrt_b(data, param[1], param[2], 0.5)
        )
    }
    optim(startval, objective, gradient, method = "L-BFGS-B", lower = c(0.01, 0.01))
}


# (Debug function)
.add_tangent_line <- function(a, b, data) {
    y <- llsbeta(data, a, b, 0.5)
    grad <- llsbeta_deriv_wrt_a(data, a, b, 0.5)
    i <- y - grad*a
    abline(i, grad)
}

