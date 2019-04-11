# Beta mixture model, based on https://stackoverflow.com/a/43561339 and code from
# Schröder, Christopher, and Sven Rahmann. 2017. Algorithms for Molecular Biology: AMB 12 (August): 21.

#' EM - EXPECTATION STEP
#' Given data, a matrix of parameters (params) and a vector of
#' mixture component proportions (pi), compute the probability
#' that each data point came from each component (weights)
#' @importFrom "matrixStats" logSumExp
#' @export
get_weights <- function(data, params, pi) {
    ncomp <- nrow(params)
    ndata <- length(data)
    stopifnot(ncomp == length(pi))
    stopifnot(all.equal(sum(pi), 1))
    ll <- matrix(0, nrow = ndata, ncol = ncomp)
    for (i in 1:ncomp) {
        ll[,i] <- log(pi[i]) + beta_ll(data, params[i,1], params[i,2])
    }
    weights <- exp(ll - apply(ll, 1, logSumExp))
    bad <- apply(weights, 1, function(row) any(is.nan(row)|is.infinite(row)))
    weights[bad & data < 0.5, ] <- rep(c(1, rep(0, ncomp-1)), each = sum(bad & data < 0.5))
    weights[bad & data > 0.5, ] <- rep(c(rep(0, ncomp-1), 1), each = sum(bad & data > 0.5))
    list(weights = weights, ll = sum(apply(ll, 1, logSumExp)))
}

#' EM - MAXIMISATION STEP
#' Given data and a matrix of mixture component weights for each
#' data point, estimate new parameter values (params) and mixture
#' proportions (pi)
get_updated_params <- function(data, weights) {
    means <- data %*% weights / colSums(weights)
    ncomp <- ncol(weights)
    vars <- sapply(1:ncomp, function(i) (data - means[i])^2 %*% weights[,i] / sum(weights[,i]))
    newpi <- colMeans(weights)
    newparams <- matrix(0, nrow = ncol(weights), 2)
    for (i in 1:ncomp) {
        newparams[i, ] <- beta_ab(means[i], vars[i])
    }
    list(pi = newpi, params = newparams)
}

#' Return the largest scaled absolute difference among all parameter
#' values in a comparison of old and new values
get_delta <- function(new_params, old_params, new_pi, old_pi) {
    epi <- max_relerror(old_pi, new_pi)
    ea <- max_relerror(old_params[, 1], new_params[, 1])
    eb <- max_relerror(old_params[, 2], new_params[, 2])
    max(epi, ea, eb)
}

#' Return the maximum scaled absolute difference between vectors
#' va and vb
max_relerror <- function(va, vb) {
    if (all(va == vb)) return (0)
    numer <- abs(va - vb)
    denom <- apply(matrix(c(va, vb), ncol = 2), 1, max)
    max(numer / denom)
}

#' Plot a single iteration of the mixture model estimation procedure of `estimate_mixture`
#' @importFrom "grid" popViewport pushViewport grid.draw
#' @importFrom "gridExtra" tableGrob ttheme_minimal
#' @importFrom "gridBase" baseViewports
#' @param "data" The observations being fitted
#' @param "params" Current Beta component `a` and `b` parameters
#' @param "pi" Component mixture weights
#' @param "step" Current iteration number
plot_estimation_step <- function(data, params, pi, step) {
    # Calculate likelihood of data given parameters
    weights_ll <- get_weights(data, params, pi)
    weights <- weights_ll$weights
    ll <- weights_ll$ll

    # Calculate data for line of uncertainty, and the boundaries
    x <- seq(0, 1, length.out = 1001)
    w <- get_weights(x, params, pi)$weights
    uncertainty <- 1 - apply(w, 1, max)

    bounds <- c()

    for (i in 2:nrow(params)) {
        bounds <- append(
            bounds,
            beta_boundary(params[i-1, 1], params[i-1, 2], params[i, 1], params[i, 2], pi[i-1], pi[i], log = FALSE)
        )
    }

    # bounds <- c(
    #     beta_boundary(params[1, 1], params[1, 2], params[2, 1], params[2, 2], pi[1], pi[2], log = FALSE),
    #     beta_boundary(params[2, 1], params[2, 2], params[3, 1], params[3, 2], pi[2], pi[3], log = FALSE)
    # )

    # Draw the histogram
    layout(matrix(c(1,1,1,3,
                    1,1,1,1,
                    2,2,2,2), nrow = 3, byrow = TRUE))
    par(mar = c(-0.1, 4, 4, 2) + 0.1)
    hist(data,
         breaks = seq(0, 1, length.out = 51),
         border = "white",
         col = hcl(240, 50, 80),
         probability = TRUE,
         main = paste("step", step, "ll =", round(ll, 4)),
         #sub = paste(round(pi, 2), collapse = " "),
         xaxt = "n", xlab = NA)

    # Add lines for the component densities
    dens <- t(t(apply(params, 1, function(r) dbeta(x, r[1], r[2]))) * pi)
    apply(dens, 2, function(c) lines(x, c, lty = 3, col = hcl(240, 80, 50)))
    lines(x, rowSums(dens), col = hcl(240, 100, 30), lwd = 2, lty = 1)

    # Add boundaries
    abline(v = bounds, col = "red")

    # Plot uncertainty
    par(mar = c(5, 4, -0.1, 2) + 0.1)
    plot(x, uncertainty, ylim = c(0, 1), type = "l", bty = "n", xlab = "VAF", ylab = "Uncertainty")
    abline(v = bounds, col = "red")

    # Draw table of parameters
    df <- data.frame(a = params[, 1], b = params[, 2], pi = pi)
    frame()
    # Grid regions of current base plot (ie from frame)
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)
    mytheme <- ttheme_minimal(
        core = list(fg_params=list(cex = 0.7)),
        colhead = list(fg_params=list(cex = 0.6)),
        rowhead = list(fg_params=list(cex = 0.6)))
    grob <- tableGrob(round(df, 2), theme = mytheme)
    grid.draw(grob)
    popViewport(3)
}

#' Run Expectation Maximisation to optimise the fit of a mixture of Beta components
#' to data
#' @param data (numeric matrix) - Observations in the range [0,1]
#' @param init_params (numeric matrix) - Matrix of initial Beta(a, b) parameters, with
#' one row per mixture component.
#' @param init_pi (numeric vector) - Vector of initial mixture component proportions,
#' with one entry per mixture component.
#' @param max_steps (integer) - Maximum number of iterations to perform while optimising
#' @param tol (numeric) - Terminate optimisation when change in parameters between iterations
#' is less than tol
#' @param plot (bool) - Plot optimisation progress to active graphical device
#' @return list(params, pi, usedsteps, assignment, uncertainty) - Optimised parameters and mixture proportions,
#' number of steps used, the assignment of each point to a mixture component,
#' and the degree of uncertainty for each point
#' @export
estimate_mixture <- function(data, init_params, init_pi,
                             max_steps = 1000, tol = 1e-5, plot = FALSE) {
    n <- length(data)
    ncomp <- nrow(init_params)

    if (!all.equal.numeric(sum(init_pi), 1)) stop(paste0("init_pi must sum to 1: (", paste(init_pi, collapse=", "), ")"))
    if (any(init_params <= 0)) stop("init_params must all be > 0")
    if (!nrow(init_params) == length(init_pi)) stop("Number of rows of init_params must equal length of init_pi")

    step <- 0
    params <- init_params
    pi <- init_pi

    if (plot) {
        plot_estimation_step(data, params, pi, step)
    }

    for (i in 1:max_steps) {
        params_old <- params
        pi_old <- pi

        w <- get_weights(data, params, pi)
        weights <- w$weights

        new <- get_updated_params(data, weights)
        params <- new$params
        pi <- new$pi

        step <- step + 1
        delta <- get_delta(params, params_old, pi, pi_old)
        if (step %% 10 == 0) print(paste(step, delta))

        if (plot) {
            plot_estimation_step(data, params, pi, step)
        }
        if (delta < tol) {
            break
        } else {
            params <- new$params
            pi <- new$pi
        }
    }
    bic <- log(length(data)) * nrow(params) - 2 * w$ll

    # Get the means of the component distributions, and sort everything by increasing mean
    means <- apply(params, 1, function(row) row[1] / sum(row))
    sorted <- order(means)
    means <- means[sorted]
    params <- params[sorted, ]
    pi <- pi[sorted]
    weights <- weights[, sorted]

    list(params = params, pi = pi, usedsteps = step,
         means = means, classification = apply(weights, 1, which.max),
         uncertainty = 1 - apply(weights, 1, max),
         data = data, ll = w$ll, bic = bic)
}

#' Get initial parameter estimates for the mixture model
#' @param v Model data
#' @param k Number of clusters
#' @return (List) Initial estimates of params and pi
#' @export
get_init <- function(v, k = 3) {
    params <- matrix(0, ncol = 2, nrow = k)
    pi <- rep(0, k)
    bounds <- c(0, 1:k/k)
    for (i in 1:k) {
        dat <- v[v>bounds[i] & v <= bounds[i+1]]
        params[i, ] <- beta_ab(mean(dat, na.rm = TRUE), var(dat, na.rm = TRUE))
        pi[i] <- length(dat) / length(v)
    }
    pi <- pi / sum(pi)
    params[params < 0.01] <- 0.01
    list(params=params, pi=pi)
}
