transform.axis <- function(data) {
# transform 2D data via their eigenvector
    e <- eigen(cov(data))

    theta <- atan(e$vectors[,2][2] / e$vectors[,2][1])
    #if (abs(theta) > abs(atan(e$vectors[,1][2] / e$vectors[,1][1]))) {
    #    theta = atan(e$vectors[,1][2] / e$vectors[,1][1])
    #}
    tm <- matrix( c(cos(theta), -sin(theta), sin(theta), cos(theta)), ncol=2,
                  nrow=2 )

    t(tm %*% t(data))
}

point.selection <- function(depth, data, alpha, lb=TRUE) {
# compare given indexes against transformed data at given alpha level
# return the true points
    if (lb) {
        which( (depth$dep < alpha) & (data[,2] < 0) == TRUE )
    } else {
        which( (depth$dep < alpha) & (data[,2] >= 0) == TRUE )
    }
}

mahalanobis.origin <- function(data, sigma) {
# calculate the mahalanobis depth with respect to the origin
    sigma.inv <- try(solve(sigma), silent = TRUE)
    if (!is.matrix(sigma.inv)) {
        sv <- svd(sigma)
        sigma.inv <- sv$v %*% diag(1/sv$d) %*% t(sv$u)
        warning("Inverse of sigma computed by SVD")
    }
    apply(data, 1, function(x) x %*% sigma.inv %*% matrix(x))
}

calculate.state_matrix <- function(nvar) {
# calculate state probabilities
    n = nrow(nvar)
    states = c()
    smatrix = c()
    for (cidx in 1:ncol(nvar)) {
        ftmp <- as.factor(nvar[,cidx])
        stmp <- levels(ftmp)
        mtmp <- matrix(0, n, length(stmp))
        sidx = 1
        for (s in stmp) {
            mtmp[which(nvar[,cidx] == s), sidx] = 1
            sidx = sidx + 1
        }
        states <- c(states, stmp)
        smatrix = cbind(smatrix, mtmp)
    }
    smatrix = as.data.frame(smatrix)
    colnames(smatrix) <- states
    smatrix
}

calculate.state_probs <- function(smatrix) {
    states = c()
    counter = c()
    for (ridx in 1:nrow(smatrix)) {
        s <- paste(smatrix[ridx,], collapse="")
        if (!(s %in% states)) {
            states = c(states, s)
            counter = c(counter, 1)
        } else {
            counter[which(states == s)] = counter[which(states == s)] + 1
        }
    }
    names(counter) <- states
    counter
}



ellipse <- function(center, sigma, n) {
# generalized equation for ellipse centered at (0, 0)
    sigma.inv <- try(solve(sigma), silent = TRUE)
    if (!is.matrix(sigma.inv)) {
        sv <- svd(sigma)
        sigma.inv <- sv$v %*% diag(1/sv$d) %*% t(sv$u)
        warning("Inverse of sigma computed by SVD")
    }


    function(mu1, mu2) {
        n * sigma.inv[1, 1] * (center[1] - mu1)^2 + n * sigma.inv[2, 2] * (center[2] - mu2)^2 + 2 * n * sigma.inv[1, 2] * (center[1] - mu1) * (center[2] - mu2)
    }
}

level <- function(n, p, alpha) {
    p * (n -1) * qf(1 - alpha, p, n-p) / (n - p)
}
