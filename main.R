# This script can be used to initiate the analysis for alphaCore, weighted 
# k-Core, and rich club analysis for defined networks
library(dplyr)


setwd("/mnt/alphaCore")


# parameters
param.data.idx <- 3 ; # index to select dataset from below
param.alphaCore.stepsize <- 0.0005;
# valid options c("k", "s") for degress and weight
param.richclub.type = "s"
param.richclub.iterations = 1000
param.richclub.seed = 1
param.richclub.reshuffle = "links"
# valid options c("alphaCore", "kCore", "richClub") 
param.analysis <- "richClub" 


data.network.fn <- c("./data/network.citation-statistics.txt",
                     "./data/network.airport-US2010.txt",
                     "./data/network.collaboration-netscience.txt",
                     "./data/network.protein-4932.small.txt")

data.labels.fn <- c("./data/label.citation-statistics.txt",
                    "./data/label.airport-US2010.txt",
                    "./data/label.collaboration-netscience.txt",
                    NULL)

data.network <- read.csv(file=data.network.fn[param.data.idx],
                         sep=" ", header=F)
data.network <- distinct(data.network)
duplicates <- which(data.network[,1] == data.network[,2])
if (!(length(duplicates) == 0)) {
  data.network = data.network[-duplicates,]
}
colnames(data.network) <- c("from", "to", "weight")

data.labels <- NULL

# load network vertex labels
if (!(is.null(data.labels.fn[param.data.idx]))) {
    data.labels <- read.csv(file=data.labels.fn[param.data.idx], header=F)
    data.labels <- as.character(data.labels$V1)
    data.labels = cbind(seq.int(data.labels), data.labels)
    data.labels = as.data.frame(data.labels)
    colnames(data.labels) <- c("id", "data.labels")
}

load.graph <- function(d) {

    g <- graph_from_edgelist(as.matrix(d[,1:2]), 
                                      directed=T)
    g <- graph_from_edgelist(as.matrix(d[,1:2]),
                                      directed=T)%>%
         set_vertex_attr("idx", value=as.character(seq(1, vcount(g))))
    E(g)$weight <- data.network$weight
    g <- simplify(g)
    g<-g%>%
      set_vertex_attr("idx", value = as.character(seq(1,vcount(g))))
    
    return(g)
}


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
      
      inputGr <- delete_vertices(inputGr, nodes_to_delete)
      nodes_to_delete <- c()
    }
  }
  cores = cbind(cores[,1], cores)
  cores = as.data.frame(cores)
  colnames(cores) <- c("node", "id", "k")
  return (cores)
}


if (param.analysis == "kCore") {
    #source("kcore.R")
    data.graph <- load.graph(data.network)
    result <- weightedkCore(data.graph)
    if (!(is.null(data.labels))) {
        result <- merge(result, data.labels, all.x=T)
    } 
} else if (param.analysis == "richClub") {
    library(tnet)
    source("richClub.R")
    rc.out <- weighted_richclub_w(data.network,
                                  rich=param.richclub.type,
                                  reshuffle=param.richclub.reshuffle, 
                                  NR=param.richclub.iterations, 
                                  seed=param.richclub.seed)
    
    
    plot(rc.out[,c("x","y")], type="b", log="x", xlim=c(1, max(rc.out[,"x"]+1)), ylim=c(0,max(rc.out[,"y"])+0.5), xaxs = "i", yaxs = "i", ylab="Weighted Rich-club Effect", xlab="Prominence (degree greater than)")
    lines(rc.out[,"x"], rep(1, nrow(rc.out)))
    data.graph <- load.graph(data.network)
    # delete vertices that iGraph fills in, or any vertices that have 0 degrees 
    # i.e. not connected to the rest of the network
    data.graph<-delete_vertices((data.graph), degree(data.graph)==0)
    
    #Store initial neighbor and weights as a record.
    initialNodeFeatures<-matrix(0, vcount(data.graph), 3)
    
    # should be faster to iterate though the dataset
    counter = 1
    for(v in as.integer(vertex_attr(data.graph)$idx)) {
      nidx <- which(data.network[,2] == v)
      inDegree <- length(nidx)
      inWeight <- sum(data.network[nidx, 3])
      initialNodeFeatures[counter,] = c(v,inDegree, inWeight) 
      counter = counter + 1
    }
    
    bins <- degree_w(rc.out)
    result <- matrix(0, vcount(data.graph), 3)
    bidx = 0
    for (b in rc.out[,1]){ 
      bidx = bidx + 1
      idx = which(initialNodeFeatures[,3] > b)
      result[idx,2:3] = cbind(initialNodeFeatures[idx,1], rep(b)) 
    }
    result[,1] = 1:vcount(data.graph)
    
    colnames(result)<-c("Id","Node","bin")
    write.csv(result, file = "richClub.csv")
}



