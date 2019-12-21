# This script can be used to initiate the analysis for alphaCore, weighted 
# k-Core, and rich club analysis for defined networks
library(dplyr)

setwd("/mnt/alphaCore")
source("kcore.R")

# parameters
param.data.idx <- 1 ; # index to select dataset from below
param.alphaCore.stepsize <- 0.0005;
# valid options c("k", "s") for degress and weight
param.richclub.type = "s"
param.richclub.iterations = 100
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
    return(g)
}


if (param.analysis == "kCore") {
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
        
}
