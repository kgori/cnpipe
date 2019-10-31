#' Implementation of PCF algorithm from 
#' Copynumber: Efficient algorithms for
#' single and multi-track copy number 
#' segmentation
#' (not very efficient version)
#' example:
#' data = c(rnorm(495, 0, 0.1), rnorm(10, 0.5, 0.1), rnorm(495, 0, 0.1), rnorm(100, 1, 0.1), rnorm(1000, 1, 0.1))
#' pcf(data)
pcf <- function(y, gamma = 40) {
    p = length(y)
    a0 = numeric(0)
    e0 = 0
    aprev = a0
    eprev = e0
    gamma = gamma * var(y)
    t <- rep(0, p) # backtrack
    for (k in seq(1, p)) {
        anow <- c(aprev, 0) + y[k]
        dnow <- anow^2 / (1:k - k - 1)
        score <- (dnow + eprev + gamma)
        min_index <- which.min(score)
        enow <- c(eprev, score[min_index])
        t[k] <- min_index
        aprev <- anow
        eprev <- enow
    }
    
    segstarts <- unique(t)
    segends <- c(segstarts[2:length(segstarts)]-1, p)
    data.frame(startpos=segstarts, endpos=segends)
}

# ys is a matrix of datapoints (rows) x samples (cols)
# y <- cbind(c(rnorm(1000, 0, 0.1), rnorm(100, 1, 0.1), rnorm(1000, 0, 0.1)), c(rnorm(1000, 0, 0.1), rnorm(100, 1, 0.1), rnorm(400, 0, 0.1), rnorm(100, 1, 0.1), rnorm(500, 0, 0.1)))
multipcf <- function(ys, gamma = 40) {
    p = nrow(ys)
    a0 = matrix(0, nrow = nrow(ys), ncol = ncol(ys))
    e0 = 0
    aprev = a0
    aprev[1,] <- ys[1,]
    anow = a0
    eprev = e0
    t <- rep(0, p) # backtrack
    for (k in seq(1, p)) {
        anow[1:(k+1)] <- aprev[1:(k+1), ] + ys[k,]
    }
}