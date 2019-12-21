# k-core decomposition given graph G
library(igraph)

### unweighted k core decomposition ###

kCore <- function(inputGr) {
    k <- 1
    nodes_to_delete <- c()
    cores <- c()

    while (vcount(inputGr) > 0) {
        for (i in V(inputGr)) {
            d <- degree(inputGr, i, mode="all")
            if (d < k) {
                nodes_to_delete <- c(nodes_to_delete, i)
            }
        }
        if (length(nodes_to_delete) == 0) {
            k <- k + 1
        } else {
            cores <- rbind(cores,
                           cbind(vertex_attr(inputGr)$idx[nodes_to_delete], k))
            inputGr <- delete.vertices(inputGr, nodes_to_delete)
            nodes_to_delete <- c()
        }
    }
    cores = cbind(cores[,1], cores)
    cores = as.data.frame(cores)
    colnames(cores) <- c("node", "id", "k")
    return (cores)
}

### weigthed k core decomposition according to Antonios Garas 2012 ###

weightedkCore <- function(inputGr) {
    k <- 1
    alpha <- 1
    beta <- 1
    nodes_to_delete <- c()
    cores <- c()

    while (vcount(inputGr) > 0) {
        for (i in V(inputGr)) {
            # get adjacent edges and sum the weights
            adj_v <- unlist(adjacent_vertices(inputGr, i, mode="in"))
            wi <- 0
            di <- length(adj_v)
            for (ai in adj_v) {
                wi <- wi + edge_attr(inputGr, "weight", 
                                     get.edge.ids(inputGr, c(ai, i)))
            }
            ki <- ((di**alpha) * (wi**beta))**(1/(alpha+beta))
            if (k > ki) {
                nodes_to_delete <- c(nodes_to_delete, i)
            }
        }
        if (length(nodes_to_delete) == 0) {
            k <- k + 1
        } else {
            message("Number of nodes to be deleted: ", length(nodes_to_delete))
            cores <- rbind(cores, 
                           cbind(vertex_attr(inputGr)$idx[nodes_to_delete], k))

            inputGr <- delete.vertices(inputGr, nodes_to_delete)
            nodes_to_delete <- c()
        }
    }
    cores = cbind(cores[,1], cores)
    cores = as.data.frame(cores)
    colnames(cores) <- c("node", "id", "k")
    return (cores)
}

### Test Graph + unweighted k-core ###
test <- function() {

    G <- graph(c(1,2,  1,3,  2,3,  2,6,  3,4,  3,5,
                 3,6,  3,7,  4,5,  4,7,  4,8,  5,7,
                 5,8,  6,7,  6,9,  7,8,  7,9), 
               directed=F)

    for (i in 1:ecount(G)) {
        G <- set_edge_attr(G, "weight", i, i)
    }

    result <- kCore(G)
}
