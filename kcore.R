# k-core decomposition given graph G
library(igraph)

### Test Graph + unweighted k-core ###
.f <- function() {

G <- graph(c(1,2,  1,3,  2,3,  2,6,  3,4,  3,5,
             3,6,  3,7,  4,5,  4,7,  4,8,  5,7,
             5,8,  6,7,  6,9,  7,8,  7,9), 
          directed=F)

for (i in 1:ecount(G)) {
    G <- set_edge_attr(G, "weight", i, i)
}

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

}

### weigthed k core decomposition according to Antonios Garas 2012 ###
data<-read.csv(file="./data/networkcitation.txt", sep="\t", header=F)
colnames(data)<-c("from","to","time","weight")
tokenGr<-graph_from_edgelist(as.matrix(data[,1:2]),directed=TRUE) 
tokenGr<-graph_from_edgelist(as.matrix(data[,1:2]),directed=TRUE)%>%
    set_vertex_attr("idx", value = as.character(seq(1,vcount(tokenGr))))
E(tokenGr)$weight<-data$weight
tokenGr<-delete_vertices((tokenGr), degree(tokenGr)==0)


k <- 1

alpha <- 1
beta <- 1

Gcopy <- tokenGr
nodes_to_delete <- c()
cores <- c()

while (vcount(Gcopy) > 0) {
    for (i in V(Gcopy)) {
        # get adjacent edges and sum the weights
        adj_v <- unlist(adjacent_vertices(Gcopy, i, mode="in"))
        wi <- 0
        di <- length(adj_v)
        for (ai in adj_v) {
            wi <- wi + edge_attr(Gcopy, "weight", get.edge.ids(Gcopy, c(ai, i)))
        }
        ki <- (di**alpha) * (wi**beta)**(1/(alpha+beta))
        if (k > ki) {
            nodes_to_delete <- c(nodes_to_delete, i)
        }
    }
    if (length(nodes_to_delete) == 0) {
        #if (vcount(Gcopy) < 200) {
        #    plot(Gcopy)
        #}
        k <- k + 1
    } else {
        message("Number of nodes to be deleted: ", length(nodes_to_delete))
        cores <- rbind(cores, cbind(vertex_attr(Gcopy)$idx[nodes_to_delete], k))

        Gcopy <- delete.vertices(Gcopy, nodes_to_delete)
        nodes_to_delete <- c()
    }
}

write.csv(cores, "results.csv")

