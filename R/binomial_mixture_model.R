
# Binomial coefficient = n!/k!(n-k)! = Γ(n+1) / Γ(k+1)Γ(n-k+1)
# log Binomial coeff = lnΓ(n+1) - lnΓ(k+1) - lnΓ(n-k+1)

#' Log bmm_likelihood of data under binomial distribution
#' @param param - Parameter of binomial distribution
#' @param data - Matrix of successes, failures in columns
#' @param binom_log_constant - Logarithm of the binomial constant for data
bmm_loglik <- function(param, data, binom_log_constant) {
    lp <- log(param)
    lq <- log(1 - param)
    data[, 1] * lp + data[, 2] * lq + binom_log_constant
}

#' Convert log-bmm_likelihood matrix to probabilities
#' @param llmat - Matrix with log-bmm_likelihoods for each mixture component in columns
#' @importFrom "Rfast" rowMaxs rowsums
bmm_ll_to_prob <- function(llmat) {
    maxs <- rowMaxs(llmat, value = TRUE)
    if (ncol(llmat) == 1) {
        probs <- matrix(1, nrow(llmat))
    } else {
        probs <- exp(llmat - maxs)
    }
    scale <- rowsums(probs)
    probs / scale
}

#' Overall log bmm_likelihood of the data given the mixture model
#' @param m Matrix of component log bmm_likelihoods
#' @importFrom "Rfast" rowMaxs rowsums
bmm_likelihood <- function(m) {
    maxs <- rowMaxs(m, value = TRUE)
    sum(log(rowsums(exp(m - maxs))) + maxs)
}

#' Fits a binomial mixture model to data in y. The number of components
#' is inferred from the length of the theta starting parameters vector.
#' @param y Data to fit, as matrix of successes and failures in columns
#' (in DNA variant context, a 'success' is an ALT read, so they go in
#' column 1)
#' @param theta Vector of starting binomial probability parameters,
#' one for each mixture component
#' @param lambda Vector of mixture weights, one for each component
#' @param tol Optimisation will stop when the bmm_likelihood improves by
#' less than this tolerance
#' @param max_iter Optimisation stops when this number of iterations is
#' reached
#' @importFrom "Rfast" rowsums Lgamma eachrow colmeans
#' @importFrom "logging" loginfo logwarn
#' @export
binomial_MM <- function(y, theta, lambda, tol = 1e-3, max_iter = 1000, fix_theta = FALSE, fix_lambda = FALSE) {
    stopifnot(length(theta) == length(lambda))
    stopifnot(all(theta > 0 & theta < 1))
    stopifnot(all(lambda > 0 & lambda < 1))
    if (sum(lambda) != 1) lambda <- lambda / sum(lambda)

    best <- -Inf

    # Log binomial coefficient
    lbc <- Lgamma(rowsums(y) + 1) - rowsums(Lgamma(y + 1))

    for (i in seq_len(max_iter)) {
        ll <- eachrow(sapply(theta, bmm_loglik, data = y, binom_log_constant = lbc),
                             log(lambda),
                             '+')
        ll[is.na(ll)] <- 0
        current <- bmm_likelihood(ll)

        if (i %% 10 == 1) loginfo("Iter = %d - likelihood = %.3f", i, current)

        p <- bmm_ll_to_prob(ll)

        if (abs(current - best) < tol) {
            loginfo("Stopping as bmm_likelihood increase is within tolerance")
            break
        }

        if (current < best) {
            logwarn("bmm_likelihood has worsened")
        } else {
            best <- current
        }

        if (i == max_iter) {
            loginfo("Stopping as max iterations have been reached")
            break
        }

        if (!fix_theta) {
            theta <- t(p) %*% y
            theta <- (theta / rowsums(theta))[, 1]
            theta <- pmin(1-1e-16, pmax(1e-16, theta))
        }
        if (!fix_lambda) {
            lambda <- colmeans(p) # No na.rm = TRUE, but is it necessary?
        }
    }

    o <- order(theta)
    list(theta = theta[o],
         lambda = lambda[o],
         probs = p[, o],
         likelihood = best,
         bic = (2 * length(theta) - 1) * log(nrow(y)) - 2 * best,
         iterations = i)
}

#' Fits a zero-inflated binomial mixture model to data in y. The number
#' of components is inferred from the length of the theta starting parameters
#' vector.
#' @param y Data to fit, as matrix of successes and failures in columns
#' (in DNA variant context, a 'success' is an ALT read, so they go in
#' column 1)
#' @param theta Vector of starting binomial probability parameters,
#' one for each mixture component. The first element should be 0, to handle
#' the zero-inflation
#' @param lambda Vector of mixture weights, one for each component
#' @param tol Optimisation will stop when the bmm_likelihood improves by
#' less than this tolerance
#' @param max_iter Optimisation stops when this number of iterations is
#' reached
#' @importFrom "Rfast" rowsums Lgamma eachrow colmeans
#' @importFrom "logging" loginfo logwarn
#' @export
zi_binomial_MM <- function(y, theta, lambda, tol = 1e-3, max_iter = 1000) {
    stopifnot(length(theta) == length(lambda))
    stopifnot(all(theta >= 0 & theta < 1))
    stopifnot(all(lambda > 0 & lambda < 1))
    if (theta[1] != 0) {
        theta[1] <- 0
    }
    if (sum(lambda) != 1) lambda <- lambda / sum(lambda)

    best <- -Inf

    # Log binomial coefficient
    lbc <- Lgamma(rowsums(y) + 1) - rowsums(Lgamma(y + 1))
    ncat <- length(theta)

    for (i in seq_len(max_iter)) {
        ll <- eachrow(sapply(theta[2:ncat], bmm_loglik, data = y, binom_log_constant = lbc),
                      log(lambda[2:ncat]),
                      '+')
        ll[is.na(ll)] <- 0
        ll <- cbind(ifelse(y[, 1] == 0, log(lambda[1]), -Inf), ll)
        current <- bmm_likelihood(ll)

        if (i %% 10 == 1) loginfo("Iter = %d - likelihood = %.3f", i, current)

        p <- bmm_ll_to_prob(ll)

        if (abs(current - best) < tol) {
            loginfo("Stopping as bmm_likelihood increase is within tolerance")
            break
        }

        if (current < best) {
            logwarn("bmm_likelihood has worsened")
        } else {
            best <- current
        }

        if (i == max_iter) {
            loginfo("Stopping as max iterations have been reached")
            break
        }

        theta <- t(p) %*% y
        theta <- (theta / rowsums(theta))[, 1]
        theta <- pmin(1-1e-16, pmax(1e-16, theta))
        theta[1] <- 0
        lambda <- colmeans(p) # No na.rm = TRUE, but is it necessary?
    }

    o <- order(theta)
    list(theta = theta[o],
         lambda = lambda[o],
         probs = p[, o],
         likelihood = best,
         bic = (2 * length(theta) - 1) * log(nrow(y)) - 2 * best,
         iterations = i)
}
