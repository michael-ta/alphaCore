library(igraph)
library(qgraph)
library(dplyr)

# Version V.0.1
# given alphaCore result output - draw the plot of the network showing the
# lowest alpha nodes + the decomposition of the network across alpha the
# various alpha levels


##### functions ######

getVertexColors<-function(inputGr, aCM) {
    vcolor <- c()
    total <- length(unique(as.numeric(aCM[,3])))
    level <- 0
    count <- 0
    first <- TRUE
    for (alpha in sort(unique(as.numeric(aCM[,3])), decreasing=F)) {
        level = level + 1
        idx  <- which(as.numeric(aCM[,3]) == alpha)
        count = count + length(idx)
        color <- rgb(1 - (count/vcount(inputGr) * (level/total)), 
                     0,
                     count/vcount(inputGr) * (level/total) )
        tidx <- which(as.numeric(vertex_attr(inputGr, "idx")) %in% 
                      as.numeric(aCM[idx,2]))
  
        if (first) {
            color <- rgb(0, 1, 0)
            first <- FALSE
        }
        vcolor[ tidx ] = color
    }
    return (vcolor)
}

######################


args <- commandArgs(trailing=F)
script.path <- sub("--file=", "", args[grep("--file=", args)])
if (!(identical(script.path, character(0)))) {
    script.basename <- dirname(script.path)
} else {
    # for debug / rstudio purposes
    script.basename <- "/mnt/alphaCore"
}
setwd(script.basename)

args <- commandArgs(trailing=T)
network <- if (!(is.na(args[2]))) args[2] else "airport-US2010"
debug.result.fn <- paste("results", network, "5e-02", "csv", sep=".")
resultfn <- if (!(is.na(args[1]))) args[1] else debug.result.fn
step.size <- sub( paste("results", network, "" ,sep="."), "", 
                  basename(resultfn) )
step.size <- sub( ".csv", "", step.size )

opt.max.node <- 1500
data.network <- read.csv(file=paste("./data/network", network, "txt",
                                    sep="."),
                         sep=" ", header=F)

labelfn <- paste("./data/label", network, "txt", sep=".")
if (! file.exists(labelfn)) {
    data.labels <- NA
} else {
    data.labels <- read.csv(labelfn, header=F)
    data.labels <- as.character(data.labels$V1)
    data.labels <- NA
}

# pre-process the data
data.network <- distinct(data.network)
duplicates <- which(data.network[,1] == data.network[,2])
if (!(length(duplicates) == 0)) {
    data.network = data.network[-duplicates,]
}

tokenGr <- graph_from_edgelist(as.matrix(data.network[,1:2]), directed=T)
tokenGr <- tokenGr %>%
    set_vertex_attr("idx", value=as.character(seq(1,vcount(tokenGr))))

tokenGr <- delete_vertices((tokenGr), degree(tokenGr)==0)

initialNodeFeatures<-matrix(0, vcount(tokenGr), 3)
  
# should be faster to iterate though the dataset
counter = 1
for(v in as.integer(vertex_attr(tokenGr)$idx)) {
    nidx <- which(data.network[,2] == v)
    inDegree <- length(nidx)
    inWeight <- sum(data.network[nidx, 3])
    initialNodeFeatures[counter,] = c(v,inDegree, inWeight) 
    counter = counter + 1
}
colnames(initialNodeFeatures)<-c("node","inDegree","inWeight")

vertex_attr(tokenGr)$label.cex = c(rep(0.5, vcount(tokenGr)))
vertex_attr(tokenGr)$label.dist = c(rep(.5, vcount(tokenGr)))
vertex_attr(tokenGr)$label.degree = c(rep(-pi/6, vcount(tokenGr)))
vertex_attr(tokenGr)$width = c(rep(0.5, vcount(tokenGr)))
vertex_attr(tokenGr)$label.color = c(rep("black", vcount(tokenGr)))
vertex_attr(tokenGr)$size = 1 + ((initialNodeFeatures[,2] / 
                            max(initialNodeFeatures[,2])) * 3)
tokenGr <- tokenGr %>% set_edge_attr("color", value=rgb(0.7, 0.7, 0.7, 0.25))

data.result <- read.csv(file=resultfn, sep=",", header=T)
data.result = data.result[,-1]

didx <- c()
alevels <- as.data.frame(data.result)
for (i in 1:15) {
    new_didx <- which( alevels[,3] == sort(unique(alevels[,3]))[i] )
    if (length(new_didx) + length(didx) < opt.max.node ) {
        didx = c(didx, new_didx)
    } else {
        didx = c(didx, new_didx[1:(opt.max.node - length(didx))])
        break
    }
}

dnodes <- as.numeric( data.result[didx, 2] )
sidx <- which(as.numeric(vertex_attr(tokenGr)$idx) %in% dnodes)
sGr <- induced_subgraph(tokenGr, sidx)
vertex_attr(sGr)$color = getVertexColors(sGr, data.result[didx,])

e <- get.edgelist(sGr)
l <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(sGr))

if (!(is.na(data.labels))) {
    data.labels=data.labels[sort(dnodes)]
}

pdf(paste("aCore.remaining", network, step.size, "pdf", sep="."))
plot(sGr,
     layout=l,
     label.dist=0.5,
     label.cex=0.05,
     vertex.label=data.labels,
     edge.arrow.size=0.1,
     width=0.25,
     rescale=T,
     asp=0)
dev.off()

colnames(alevels) <- c("id", "node", "alpha")
alevels <- as.data.frame(alevels) %>% 
           group_by(alpha) %>% 
           summarise(no_rows = length(alpha))
oidx <- order(as.numeric(alevels$alpha), decreasing=T)
alevels.core <- cbind(cumsum(alevels$no_rows[oidx]) / sum(alevels$no_rows),
                      alevels$alpha[oidx])

pdf(paste("aCore.decompose", network, step.size, "pdf", sep="."))
plot(alevels.core[,1], alevels.core[,2], type="b", cex=1, pch=18,
     xlab="percentage of nodes", ylab="alpha level")
dev.off()

count <- 0
pdf(paste("aCore.degree_vs_weights", network, step.size, "pdf", sep="."))
plot(initialNodeFeatures[,2:3], pch=".")
first <- T
level <- 0
total <- length(alevels$alpha)

for (alpha in sort(alevels$alpha, decreasing=F)) {
    level = level + 1
    idx <- which(as.numeric(data.result[,3]) == alpha)

    count = count + length(idx)
    color <- rgb(1 - (count/vcount(tokenGr) * (level/total)),
                 0,
                 count/vcount(tokenGr) * (level/total))
    tidx <- which(as.numeric(vertex_attr(tokenGr, "idx")) %in% 
                  as.numeric(data.result[idx,2]))
    if (first) {
        color <- rgb(0, 1, 0)
        first <- F
    }
    points(initialNodeFeatures[tidx, 2:3], col=color, pch=".")    
}
dev.off()
