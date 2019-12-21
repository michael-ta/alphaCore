library(tnet)

richclubPlot <- function(out) {
    plot(out[,c("x", "y")], type="b", log="x", 
         xlim=c(1, max(out[,"x"])+1),
         ylim=c(0, max(out[,"y"])+0.5), xaxs="i",
         ylab="weighted Rich-Club Effect", 
         xlab="Prominence (weight greater than)")
    lines(out[, "x"], rep(1, nrow(out))) 
}

richclubCore <- function(d, levels=NULL) {
    if (is.null(levels)) {
        stop("invalid level provided")
    }

    cores <- c()
    # currently we focus on the weights
    for (l in sort(levels, decreasing=T)) {
        rc.nodes <- which(aggregate( . ~ from, data=d, FUN=sum)[,3] > l)
        vidx <- aggregate( . ~ from, data=d, FUN=sum)[rc.nodes,]$from
        
        if (!(length(cores) == 0)) {
            vidx <- setdiff(vidx, cores[,1])
        }
        if (length(vidx) > 0) {
            cores = rbind(cores,
                          cbind(vidx, l) )
        }
    }

    cores = cbind(seq(1, length(cores[,1])), cores)
    cores = as.data.frame(cores)
    colnames(cores) <- c("node", "id", "k")

    return (cores)
}

