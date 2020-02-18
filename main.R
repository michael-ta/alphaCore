# This script can be used to initiate the analysis for alphaCore, weighted 
# k-Core, and rich club analysis for defined networks
library(dplyr)
library(igraph)

#setwd("D:/repos/alphaCore")

# command line params
args <- commandArgs(trailing=F)

# parameters / set value for debug purposes
# valid options c("alphaCore", "kCore", "richClub") 
param.analysis <- if (!is.na(args[2])) args[2] else "kCore"
param.data.network <- if (!is.na(args[3])) args[3] else "BAT"
#param.alphaCore.stepsize <- 0.0005;
# valid options c("k", "s") for degress and weight
param.richclub.type = "s"
param.richclub.iterations = 1000
param.richclub.seed = 1
param.richclub.reshuffle = "links"

data.network.fn <- paste("./data/network.",
                         param.data.network, ".txt", sep="")
data.labels.fn <- paste("./data/label.", 
                        param.data.network, ".txt", sep="")

data.network <- read.csv(file=data.network.fn,
                         sep=" ", header=F)
data.network <- distinct(data.network)
duplicates <- which(data.network[,1] == data.network[,2])
if (!(length(duplicates) == 0)) {
  data.network = data.network[-duplicates,]
}
colnames(data.network) <- c("from", "to", "weight")

data.labels <- NA

# load network vertex labels
if (file.exists(data.labels.fn)) {
    data.labels <- read.csv(file=data.labels.fn, header=F)
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
  min_k <- NA
  
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
      # normalize node weight
      if (di > 0) {
          wi <- wi / di
      } else {
          wi <- 0
      }
      ki <- ((di**alpha) * (wi**beta))**(1/(alpha+beta))
      if (is.nan(ki)) {
        browser()
      }
      if (k > ki) {
        nodes_to_delete <- c(nodes_to_delete, i)
      }
      if (is.na(min_k)) {
        min_k = ki 
      } else if (ki < min_k) {
        min_k = ki
      }
    }
    if (length(nodes_to_delete) == 0) {
      if (min_k > k) {
          k = floor(min_k) + 1
      } else {
          k = k + 1
      }
      message(k)
      min_k = NA
    } else {
      message("Number of nodes to be deleted: ", length(nodes_to_delete),
              " for k=", k)
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
    if (!(is.na(data.labels))) {
        result <- merge(result, data.labels, all.x=T)
    } 
    colnames(result)<-c("Id","Node","bin")
    write.table(result, file = "wkCore.csv", row.names=NA)
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
    write.table(result, file = "richClub.csv", row.names=NA)
}



