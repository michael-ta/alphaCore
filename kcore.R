# k-core decomposition given graph G
library(igraph)


G <- graph(c(1,2,  1,3,  2,3,  2,6,  3,4,  3,5,
             3,6,  3,7,  4,5,  4,7,  4,8,  5,7,
             5,8,  6,7,  6,9,  7,8,  7,9), 
          directed=F)

### unweighted k core decomposition ###
k <- 1
nodes_to_delete <- c()
Gcopy <- G
while (vcount(Gcopy) > 0) {
    for (i in V(Gcopy)) {
        d <- degree(Gcopy, i, mode="all")
        if (d < k) {
            nodes_to_delete <- c(nodes_to_delete, i)
        }
    }
    if (length(nodes_to_delete) == 0) {
        plot(Gcopy)
        k <- k + 1
    } else {
        Gcopy <- delete.vertices(Gcopy, nodes_to_delete)
        nodes_to_delete <- c()
    }
}

### weigthed k core decomposition according to Antonios Garas 2012 ###
k <- 1
Gcopy <- G
for (i in 1:ecount(G)) {
    G <- set_edge_attr(G, "weight", i, i)
    Gcopy <- set_edge_attr(G, "weight", i, i)
}
nodes_to_delete <- c()
while (vcount(Gcopy) > 0) {
    for (i in V(Gcopy)) {
        # get adjacent edges and sum the weights
        adj_v <- unlist(adjacent_vertices(Gcopy, i, mode="in"))
        wi <- 0
        di <- length(adj_v)
        for (ai in adj_v) {
            wi <- wi + edge_attr(Gcopy, "weight", get.edge.ids(Gcopy, c(i, ai)))
        }
        ki <- sqrt(di * wi)
        if (ki < k) {
            nodes_to_delete <- c(nodes_to_delete, i)
        }
    }
    if (length(nodes_to_delete) == 0) {
        plot(Gcopy)
        k <- k + 1
    } else {
        Gcopy <- delete.vertices(Gcopy, nodes_to_delete)
        nodes_to_delete <- c()
    }
}


