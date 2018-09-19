# Helper functions for working with the Beta distribution

#' Find boundary between two Beta distributions
beta_boundary <- function(a1, b1, a2, b2, p1, p2, log = TRUE, plot = FALSE) {
    # Look at region between the modes of the two Beta distributions.
    # If mode finding fails, fall back to distribution means.
    mode1 <- tryCatch(beta_mode(a1, b1), error = function(e) {
        warning(paste("Fall back to Beta mean for parameters a =", a1, "b =", b1))
        beta_mean(a1, b1)
    })
    mode2 <- tryCatch(beta_mode(a2, b2), error = function(e) {
        warning(paste("Fall back to Beta mean for parameters a =", a2, "b =", b2))
        beta_mean(a2, b2)
    })
    psum <- p1 + p2
    p1 <- p1 / psum
    p2 <- p2 / psum
    x <- seq(mode1, mode2, length.out = 1001)
    if (log) {
        df <- data.frame(dens1 = log(p1) + dbeta(x, a1, b1, log=log),
                         dens2 = log(p2) + dbeta(x, a2, b2, log=log))
        y <- apply(df, 1, function(row) {
            log(abs(diff(exp(row - max(row))))) + max(row) # Log-AbsDiff-Exp function
        })
    } else {
        df <- data.frame(dens1 = p1 * dbeta(x, a1, b1, log=log),
                         dens2 = p2 * dbeta(x, a2, b2, log=log))
        y <- apply(df, 1, function(row) abs(diff(row)))
    }
    result <- x[which.min(y)]
    if (plot) {
        plot(x, y, type = "l", ylim = c(min(c(y, df$dens1, df$dens2)), max(y)))
        lines(x, df$dens1, lty=2)
        lines(x, df$dens2, lty=2)
        abline(v = result, col = "red")
    }
    result
}

#' Find mode of beta distribution
beta_mode <- function(a, b) {
    if (a > 1 & b > 1) {
        return((a - 1) / (a + b - 2))
    }
    else if (a == 1 & b == 1) result <- 0.5 # Any value in [0, 1]
    else if (a == 1 & b > 1) result <- 0
    else if (a > 1 & b == 1) result <- 1
    else stop(paste("Beta distribution has no mode: a =", a, "b =", b))
    result
}

#' Find mean of beta distribution
beta_mean <- function(a, b) {
    if (a < 0 | b < 0) stop(paste("Beta parameters out of range: a =", a, "b =", b))
    a / (a + b)
}

#' Calculate log-likelihood of data given Beta(a, b)
beta_ll <- function(data, a, b) {
    dbeta(data, a, b, log = TRUE)
}

#' Find a and b s.t. Mean(Beta(a,b)) = mn and Var(Beta(a,b)) = va
beta_ab <- function(mn, va) {
    a <- ((1-mn)/va - 1/mn) * mn^2
    b <- a * (1/mn - 1)
    c(a, b)
}
