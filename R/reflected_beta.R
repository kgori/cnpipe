# Density, density-derivative, log-likelihood and derivative-log-likelihood functions
# for a Beta distribution reflected at 0.5

#' Reflected beta distribution density function
#' @export
dreflbeta <- function(x, a, b, log = FALSE) {
    d1 <- dbeta(x, a, b, log = log)
    d2 <- dbeta(x, b, a, log = log)
    if (log) {
        densm <- matrix(c(d1, d2), ncol = 2)
        return(apply(densm, 1, function(r) {
            m <- max(r)
            log(sum(exp(r - m))) + m
        }))
    }
    d1 + d2
}

dreflbeta_deriv_wrt_x <- function(x, a, b) {
    numer1 <- x^(a - 1) * (1 - b) * (1 - x)^(b - 1)
    numer2 <- x^(a - 1) * (1 - a) * (1 - x)^(b - 1)
    numer3 <- x^(b - 1) * (1 - a) * (1 - x)^(a - 1)
    numer4 <- x^(b - 1) * (1 - b) * (1 - x)^(a - 1)
    denom1 <- (1 - x) * beta(a, b)
    denom2 <- x * beta(a, b)
    numer1 / denom1 - numer2 / denom2 + numer3 / denom1 - numer4 / denom2
}

dreflbeta_deriv_wrt_a <- function(x, a, b) {
    elem1 <- x^(a - 1) * (1 - x)^(b - 1) * (log(x) + digamma(a + b) - digamma(a))
    elem2 <- x^(b - 1) * (1 - x)^(a - 1) * (log(1 - x) + digamma(a + b) - digamma(a))
    return((elem1 + elem2) / beta(a, b))
}

dreflbeta_deriv_wrt_b <- function(x, a, b) {
    elem1 <- x^(a - 1) * (1 - x)^(b - 1) * (log(1 - x) - digamma(b) + digamma(a + b))
    elem2 <- x^(b - 1) * (1 - x)^(a - 1) * (log(x) - digamma(b) + digamma(a + b))
    return((elem1 + elem2) / beta(a, b))
}

#' Log-likelihood of a sample conditioned on these parameters, assuming reflected beta model
#' @export
ll_reflbeta <- function(params, data) {
    sum(dreflbeta(data, params[1], params[2], log = TRUE))
}

dll_reflbeta <- function(params, data) {
    ll <- dreflbeta(data, params[1], params[2], log = TRUE)
    mll <- max(ll)
    dll_wrt_a <- dreflbeta_deriv_wrt_a(data, params[1], params[2])
    dll_wrt_b <- dreflbeta_deriv_wrt_b(data, params[1], params[2])

    c(sum(dll_wrt_a / exp(ll - mll) / exp(mll)),
      sum(dll_wrt_b / exp(ll - mll) / exp(mll)))
}

#' Estimate parameters of reflected Beta distribution by maximum likelihood
#' @importFrom "nloptr" slsqp
#' @export
est_reflbeta <- function(data, startval = c(1, 10)) {
    objective <- function(param) {
        -ll_reflbeta(param, data)
    }

    gradient <- function(param) {
        -dll_reflbeta(param, data)
    }

    ineq <- function(x) {
        x[2] - x[1]
    }

    if (ineq(startval) < 0) {
        startval <- rev(startval)
    }

    slsqp(startval, objective, gradient, hin = ineq, lower = c(0.01, 0.01))
}

# (Debug function)
.add_tangent_line <- function(a, b, data, index) {
    y <- ll_reflbeta(c(a, b), data)
    grad <- dll_reflbeta(c(a, b), data)[index]
    i <- y - grad*x
    abline(i, grad)
}
